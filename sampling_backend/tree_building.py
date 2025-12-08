import json, csv
import time
from classes_ranks_definition import * # imports TaxonNode, GenomeRecord and everything needed for the rank
from enum import Enum
from collections import deque # double-ended queue, can remove elements from both sides of a list
from collections import defaultdict # provides a default value for a nonexistent key in the dictionary


'''---tree construction---'''
def parse_taxonomy(taxonomy_file): # reads .jsonl file containing taxonomy data and builds tree of TaxonNode objects
    taxon_nodes = {} # dictionary holding each node keyed by tax_id
    with open(taxonomy_file, "r") as f:
        for line in f:
            if not line.strip():
                continue  # skip empty lines
            node = json.loads(line) # each line is one taxonomic record
            taxonomy = node.get("taxonomy", {})
            tax_id = taxonomy.get("tax_id")
            if not tax_id:
                continue  # skip invalid entries
            parent_ids = taxonomy.get("parents", [])
            if parent_ids and parent_ids[-1] != tax_id: # the last parent ID in the list is the immediate parent
                parent_id = parent_ids[-1]
            else:
                parent_id = 1 # if it does not have a parent assign it to root node
            rank_str = taxonomy.get("rank")

            # normalize SUPERKINGDOM and REALM as DOMAIN
            if rank_str in {"SUPERKINGDOM", "REALM"}:
                rank = RankType.DOMAIN
                original_rank = rank_str
            else: # converts string rank into RankType Enum, if not possible NO_RANK
                original_rank = rank_str
                try:
                    rank = RankType(rank_str)
                except (ValueError, TypeError):
                    rank = RankType.NO_RANK

            name = taxonomy.get("current_scientific_name", {}).get("name", "")
            taxon_nodes[tax_id] = TaxonNode(tax_id, parent_id, original_rank, rank, name)
            taxon_nodes[tax_id].original_name = name # remember original, uninterpolated name for user lookups

        # flag basic anomalies
        orphans, self_parents = [], []
        for tax_id, n in taxon_nodes.items():
            if tax_id == 1:
                continue
            pid = getattr(n, "parent_id", None)
            n.is_self_parent = (pid == tax_id)
            if n.is_self_parent:
                self_parents.append(tax_id)
            # parent missing from the map -> orphan
            if pid is None or pid not in taxon_nodes:
                n.is_orphan = True
                orphans.append(tax_id)
            else:
                n.is_orphan = False

    # link children to parents after all nodes are created
    for tax_id, node in taxon_nodes.items():
        if tax_id == 1:
            continue # skip root from being attached as child
        parent = taxon_nodes.get(node.parent_id)
        if parent:
            parent.add_child(node)

    # determine root node
    root_tax_id = 1  # NCBI universal root
    root_node = taxon_nodes.get(root_tax_id)

    # nodes not reachable from root are unreachable
    reachable = set()
    if root_node:
        stack = [root_node]
        while stack:
            cur = stack.pop()
            if cur.tax_id in reachable:
                continue
            reachable.add(cur.tax_id)
            stack.extend(cur.children)

    unreachable = []
    for tax_id, n in taxon_nodes.items():
        if tax_id == root_tax_id:
            n.is_unreachable = False
            continue
        n.is_unreachable = (tax_id not in reachable)
        if n.is_unreachable:
            unreachable.append(tax_id)

    # stash anomalies summary on root for downstream steps
    if not root_node: # if the file was empty or no root present, return
        return None, taxon_nodes
    root_node.anomalies_total = {
        "orphans": orphans,
        "self_parents": self_parents,
        "unreachable": unreachable,
        "n_orphans": len(orphans),
        "n_self_parents": len(self_parents),
        "n_unreachable": len(unreachable),
        "sample_orphans": orphans[:5],
        "sample_self_parents": self_parents[:5],
        "sample_unreachable": unreachable[:5]}

    return root_node, taxon_nodes

'''---genome handling---'''
def attach_genomes_to_nodes(genome_file, taxon_nodes): # reads genome records from .jsonl file and associates them with taxon nodes
    with open(genome_file, "r") as f:
        for line in f:
            genome = json.loads(line) # parses each line
            tax_id = genome.get("organism", {}).get("tax_id")
            if not tax_id or tax_id not in taxon_nodes:
                continue
            # nested fields
            assembly_level = (genome.get("assembly_info", {}).get("assembly_level")
                    or genome.get("assembly", {}).get("assembly_level")
                    or genome.get("assembly", {}).get("level"))

            # refseq_category is optional and often missing
            refseq_category = (genome.get("refseq_category")
                    or genome.get("refseq", {}).get("category"))

            source_database = genome.get("source_database")  # REFSEQ or GENBANK
            # if missing but source_database == REFSEQ, mark as "reference genome"
            if not refseq_category and source_database == "SOURCE_DATABASE_REFSEQ":
                refseq_category = "reference genome"

            is_reference = (refseq_category == "reference genome") # “reference genome” flag only if refseq_category explicitly says so

            record = GenomeRecord(
                accession=genome.get("accession"),
                tax_id=tax_id,
                name=genome.get("organism", {}).get("organism_name", ""),
                refseq_category=refseq_category,  # may be None
                assembly_level=assembly_level,
                is_reference=is_reference, # only when category present
                origin_source="NCBI",
                source_database=source_database)

            taxon_nodes[tax_id].add_genome(record)

def attach_ensembl_tsv_to_nodes(tsv_files, taxon_nodes):
    required_cols = {"taxonomy_id", "assembly_accession", "assembly", "division", "core_db"}

    ok = 0
    skipped_no_tax = 0
    skipped_no_acc = 0
    skipped_missing_taxon = 0

    for path in (tsv_files or []):
        with open(path, "r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            headers = [h.strip() for h in reader.fieldnames or []] # normalize header key for display name
            # 'name' can appear as '#name' in some divisions, so normalize to 'name'
            name_key = "name" if "name" in headers else ("#name" if "#name" in headers else None)

            # header check (warn but continue)
            missing = [c for c in required_cols if c not in headers]
            if missing:
                print(
                    f"[WARN] Ensembl TSV missing required columns {missing} in {path}; will fill as empty when absent.")

            for row in reader:
                try:
                    tax_id = int((row.get("taxonomy_id") or "").strip())
                except ValueError:
                    tax_id = None
                if not tax_id:
                    skipped_no_tax += 1
                    continue

                if tax_id not in taxon_nodes:
                    skipped_missing_taxon += 1
                    continue

                accession = (row.get("assembly_accession") or "").strip()
                if not accession:
                    skipped_no_acc += 1
                    continue

                assembly = (row.get("assembly") or "").strip()
                division = (row.get("division") or "").strip()
                core_db = (row.get("core_db") or "").strip()
                name = (row.get(name_key) or "").strip() if name_key else ""

                record = GenomeRecord(
                    accession=accession,
                    tax_id=tax_id,
                    name=name,
                    refseq_category=None,
                    assembly_level=None,
                    is_reference=False,
                    origin_source="Ensembl",
                    ensembl_division=division,
                    assembly_name=assembly,
                    ensembl_core_db=core_db)

                taxon_nodes[tax_id].add_genome(record)
                ok += 1


def collect_genomes(node, source="NCBI"): # sort genomes and record their positions in linear array
    linear_genomes = []  # list holding all genome records in traversal order
    genome_positions = {}  # dict mapping tax_id to tuple (start_index, end_index) in linear_genomes list
    index = [0]  # list used for tracking running index across recursive calls

    def dfs(n):  # preorder traversal
        if n.genomes:
            genomes = [g for g in n.genomes if getattr(g, "origin_source", "NCBI") == source] # filter genomes based on source flag
            if genomes:
                start = index[0]  # save current index = position where this node's genomes start in the list
                genomes.sort(key=lambda g:
                    (not g.is_reference,
                        g.assembly_level not in ["Complete Genome", "Chromosome"],
                        g.accession))

                linear_genomes.extend(genomes)  # adds sorted genomes to global list
                end = index[0] + len(genomes)
                genome_positions[n.tax_id] = (start, end)  # record of span of genome positions for this node
                index[0] = end  # update index to continue from new position

        for child in n.children:  # go through all children of node
            dfs(child)

    dfs(node)
    return linear_genomes, genome_positions


'''---rank interpolation and completion ---'''
def annotate_tree_with_ranks_and_spans(root, genome_positions):
    # computes genome span ranges, nearest descendant rank in subtree and rank interpolation; additionally do accession discrimination of NCBI and Ensembl
    ncbi_accs = set()
    ensembl_accs = set()

    def dfs(node): # postorder
        nonlocal ncbi_accs, ensembl_accs
        # update accession sets from this node's genomes
        genomes = getattr(node, "genomes", None) or []
        for g in genomes:
            acc = getattr(g, "accession", None)
            if not acc:
                continue
            src = getattr(g, "origin_source", "NCBI")
            if src == "NCBI":
                ncbi_accs.add(acc)
            elif src == "Ensembl":
                ensembl_accs.add(acc)

        ranked_children = []
        unranked_children = []
        starts, ends = [], []
        nearest_level = None # track nearest ranked descendant (smallest rank_level below this node)

        for child in node.children: # postorder (children first)
            dfs(child)

            if hasattr(child, "genome_span") and child.genome_span != (None, None): # child genome index ranges to build span for current node
                s, e = child.genome_span
                starts.append(s)
                ends.append(e)

            # classify child by rank
            if child.rank in RANK_ORDER:
                ranked_children.append(child)
            elif child.rank in UNRANKED_TYPES:
                unranked_children.append(child)

            # track nearest ranked descendant
            if child.rank_level is not None: # consider child's own rank_level
                nearest_level = child.rank_level if nearest_level is None else min(nearest_level, child.rank_level)
            child_nearest = getattr(child, "nearest_descendant_rank_level", None) # and child's precomputed nearest
            if child_nearest is not None:
                nearest_level = child_nearest if nearest_level is None else min(nearest_level, child_nearest)

        node.nearest_descendant_rank_level = nearest_level
        node.nearest_descendant_rank = (RANK_ORDER[nearest_level] if nearest_level is not None else None)

        # assign genome span
        if node.tax_id in genome_positions:
            s, e = genome_positions[node.tax_id]
            starts.append(s)
            ends.append(e)
        node.genome_span = (min(starts), max(ends)) if starts and ends else (None, None) # computes total range covered by node and its children

        # interpolate unranked children sitting between two known ranks
        if node.rank in RANK_ORDER and ranked_children and unranked_children:
            parent_index = RANK_ORDER.index(node.rank)
            for unranked in unranked_children:
                child_nearest = getattr(unranked, "nearest_descendant_rank_level", None)
                if child_nearest is None:
                    continue  # nothing ranked below
                next_level = parent_index + 1
                if child_nearest > next_level:
                    assigned_rank = RANK_ORDER[next_level]
                    unranked.rank = assigned_rank
                    unranked.rank_level = next_level
                    if not getattr(unranked, "original_name", None):
                        unranked.original_name = unranked.name
                    unranked.name = format_interpolated_name(assigned_rank, unranked.name)
                    unranked.interpolated = True


        # handle unranked children with genomes but no descendants
        if node.rank in RANK_ORDER:
            parent_index = RANK_ORDER.index(node.rank)
            next_level = parent_index + 1
            for unranked in unranked_children:
                if unranked.rank in RANK_ORDER: # already interpolated
                    continue
                if unranked.genomes: # only interpolate if this node has genomes
                    child_nearest = getattr(unranked, "nearest_descendant_rank_level", None)
                    if child_nearest is None or child_nearest > next_level:
                        assigned_rank = RANK_ORDER[next_level]
                        unranked.rank = assigned_rank
                        unranked.rank_level = next_level
                        if not getattr(unranked, "original_name", None):
                            unranked.original_name = unranked.name
                        unranked.name = format_interpolated_name(assigned_rank, unranked.name)
                        unranked.interpolated = True

    dfs(root)

    # after full traversal, compute shared set and stash on root
    shared_accs = ncbi_accs & ensembl_accs
    root.genome_source_stats = {
        "ncbi_accessions": ncbi_accs,
        "ensembl_accessions": ensembl_accs,
        "shared_accessions": shared_accs}


def format_interpolated_name(rank, original_name):
    abbrev = RANK_ABBREVIATIONS.get(rank, rank.name[:4])
    return f"{abbrev}_{original_name}"


def fill_missing_ranks(root, taxon_nodes): # fill missing taxonomic levels by inserting phantom nodes between parent and child when their ranks are not directly adjacent
    next_id = [10_000_000_000]
    phantom_count = 0

    def insert_phantom_chain(parent, child):
        nonlocal phantom_count
        missing_levels = child.rank_level - parent.rank_level - 1 # checks how many rank levels are missing between parent and child
        if missing_levels <= 0:
            return [child]  # if no levels missing return child as it is
        if child.phantom:
            return [child]  # do not insert more phantoms beneath a phantom

        current_parent = parent
        phantom_chain = []
        for i in range(1, missing_levels + 1): # insert one phantom node per missing rank level
            rank_level = parent.rank_level + i
            rank = RANK_ORDER[rank_level]
            new_id = next_id[0]
            next_id[0] += 1

            abbrev = RANK_ABBREVIATIONS.get(rank, rank.name[:4])
            phantom_name = f"{abbrev}_of_{child.name}_{child.tax_id}"

            phantom_node = TaxonNode(
                tax_id=new_id,
                parent_id=current_parent.tax_id,
                original_rank=None,
                rank=rank,
                name=phantom_name,
                phantom=True)
            phantom_node.rank_level = rank_level
            phantom_node.genome_span = getattr(child, "genome_span", (None, None)) # gets genome span from child
            taxon_nodes[new_id] = phantom_node # add new phantom node to global list
            if phantom_node not in current_parent.children:
                current_parent.children.append(phantom_node) # makes it a child of last node in chain
            current_parent = phantom_node
            phantom_count += 1
            phantom_chain.append(phantom_node)

        #print(f"Inserting chain from {parent.name} ({parent.rank}) to {child.name} ({child.rank}): {[p.rank.name for p in phantom_chain]}")

        if child not in current_parent.children:
            current_parent.children.append(child) # attach original child to ensure chain end at original child
        return phantom_chain + [child]

    queue = deque([root])

    while queue: # breadth‑first search: rebuild each node’s children list exactly once
        node = queue.popleft()
        original_children = list(node.children)  # take a snapshot
        node.children.clear()  # prepare to rebuild it

        for child in original_children:
            if getattr(child, "phantom", False) or node.rank_level is None or child.rank_level is None: # skip phantom nodes or invalid ranks
                node.children.append(child)
                queue.append(child)
                continue

            inserted_chain = insert_phantom_chain(node, child) # insert phantom chain between node and child
            if inserted_chain[0] not in node.children:
                node.children.append(inserted_chain[0])# add the first node in phantom chain to current node
            queue.append(inserted_chain[-1]) # only queue the real (final) child to avoid re-processing chain elements

    print(f"Inserted {phantom_count} phantom nodes.")


'''---linearize tree---'''
def linearize_tree(root): # full tree traversal, returning list of leaf nodes, list of all nodes and assigns leaf span to each node
    ordered_leaves = []
    ordered_nodes = []

    def dfs(node):
        ordered_nodes.append(node)  # collect internal nodes
        start = len(ordered_leaves) # stores current index in ordered_leaves
        if not node.children: # if no children, then it is a leaf
            ordered_leaves.append(node)
        else:
            for child in node.children:
                dfs(child)
        end = len(ordered_leaves)
        node.span = (start, end)

    dfs(root)
    return ordered_leaves, ordered_nodes

'''---siblings---'''
def link_global_rank_siblings(ordered_nodes, rank_attr="rank"): # links taxonomic nodes with same rank into sibling chain
    last_seen_by_rank = {} # dict to remember last node seen for each rank (key= rank and value is node object)
    for node in ordered_nodes:
        rank = getattr(node, rank_attr, None)
        if rank is not None:
            if rank in last_seen_by_rank: # if a previous node of same rank was seen: that previous node's sibling attribute is set to current node
                last_seen_by_rank[rank].sibling = node
            last_seen_by_rank[rank] = node

'''---save tree---'''
'''
class CustomEncoder(json.JSONEncoder): # to make python objects into valid JSON
    def default(self, o): # if it encounters an unknown object
        if isinstance(o, Enum): # if object is enum it returns the name
            return o.name
        elif hasattr(o, "__dict__"): # if object is custom class handle with serialize_object()
            return self.serialize_object(o)
        elif isinstance(o, (list, tuple)): # lists and tuples are handeled for each element calling the default
            return [self.default(item) for item in o]
        elif isinstance(o, dict): # for dicts key and value are serialized
            return {key: self.default(value) for key, value in o.items()}
        elif isinstance(o, (str, int, float, bool, type(None))): # primitive values shouldn't go to .default again
            return o
        return super().default(o) # if none of the above apply let python's encoder handle it

    def serialize_object(self, obj):
        result = {}
        for key, value in obj.__dict__.items():
            result[key] = self.default(value)
        return result


def serialize_node(node): # convert a TaxonNode to a fully serializable dictionary
    data = node.__dict__.copy()
    # handle known Enum fields explicitly
    if isinstance(data.get("rank"), Enum):
        data["rank"] = data["rank"].name
    if isinstance(data.get("original_rank"), Enum):
        data["original_rank"] = data["original_rank"].name
    if "sibling" in data and hasattr(data["sibling"], "tax_id"):
        data["sibling"] = data["sibling"].tax_id
    data["children"] = [child.tax_id for child in node.children]
    return data

def serialize_genome(genome):
    data = genome.__dict__.copy()
    if isinstance(data.get("rank"), Enum):
        data["rank"] = data["rank"].name
    return data

def save_tree_output(output_prefix, taxon_nodes, linear_genomes, ordered_nodes):
    tax_file = f"{output_prefix}_taxonomy.jsonl"
    genome_file = f"{output_prefix}_genomes.jsonl"
    ordered_file = f"{output_prefix}_ordered_nodes.jsonl"
    linear_genome_file = f"{output_prefix}_linear_genomes.jsonl"
    taxon_nodes_file = f"{output_prefix}_taxon_nodes.json"

    # save all taxa with all attributes
    with open(tax_file, "w") as f:
        for node in taxon_nodes.values():
            record = serialize_node(node)
            f.write(json.dumps(record, cls=CustomEncoder) + "\n")

    # save all genome records with all attributes
    with open(genome_file, "w") as f:
        for genome in linear_genomes:
            record = serialize_genome(genome)
            f.write(json.dumps(record, cls=CustomEncoder) + "\n")

    # save ordered_nodes as list
    with open(ordered_file, "w") as f:
        full_list = [serialize_node(node) for node in ordered_nodes]
        json.dump(full_list, f, cls=CustomEncoder, indent=2)

    # save linear_genomes as list
    with open(linear_genome_file, "w") as f:
        full_genome_list = [serialize_genome(g) for g in linear_genomes]
        json.dump(full_genome_list, f, cls=CustomEncoder, indent=2)

    # save taxon_nodes as dict {tax_id: serialized node}
    with open(taxon_nodes_file, "w") as f:
        taxon_dict = {tax_id: serialize_node(node) for tax_id, node in taxon_nodes.items()}
        json.dump(taxon_dict, f, cls=CustomEncoder, indent=2)

    print(f"[Output] Saved full taxonomy to {tax_file}")
    print(f"[Output] Saved full genomes to {genome_file}")
    print(f"[Output] Saved ordered nodes list to {ordered_file}")
    print(f"[Output] Saved linear genomes to {linear_genome_file}")
    print(f"[Output] Saved taxon_nodes dict to {taxon_nodes_file}")
'''

'''---build tree basic (without phantom nodes)---'''
def build_tree_basic(taxonomy_file, genome_file, output_prefix="output_prefix_basic", ensembl_tsv_files=None):
    start_total = time.time()

    # Step 1: Parse taxonomy
    t0 = time.time()
    root_node, taxon_nodes = parse_taxonomy(taxonomy_file)
    print(f"[Step 1] Parsed taxonomy in {time.time() - t0:.2f} seconds")

    # Step 2: Attach genomes
    t0 = time.time()
    attach_genomes_to_nodes(genome_file, taxon_nodes)

    # Attach Ensembl genomes
    if ensembl_tsv_files:
        attach_ensembl_tsv_to_nodes(ensembl_tsv_files, taxon_nodes)

    print(f"[Step 2] Attached genomes in {time.time() - t0:.2f} seconds")

    # Step 3: Collect genomes and genome positions
    t0 = time.time()
    linear_genomes, genome_positions = collect_genomes(root_node)
    print(f"[Step 3] Collected genome positions in {time.time() - t0:.2f} seconds")

    # Step 4: Interpolate ranks and assign genome spans
    t0 = time.time()
    annotate_tree_with_ranks_and_spans(root_node, genome_positions)
    print(f"[Step 4] Annotated tree (interpolation & spans) in {time.time() - t0:.2f} seconds")

    # Step 5: Linearize tree
    t0 = time.time()
    ordered_leaves, ordered_nodes = linearize_tree(root_node)
    print(f"[Step 6] Linearized tree in {time.time() - t0:.2f} seconds")

    # Step 6: Link siblings
    t0 = time.time()
    link_global_rank_siblings(ordered_nodes)
    print(f"[Step 7] Linked siblings in {time.time() - t0:.2f} seconds")

    # Step 7: Build name index
    t0 = time.time()
    name_to_taxids = defaultdict(list)
    for node in taxon_nodes.values():
        orig = getattr(node, "original_name", None)
        if orig:
            name_to_taxids[orig.lower()].append(node.tax_id)
    print(f"[Step 8] Built name index in {time.time() - t0:.2f} seconds")

    print(f"\n[Total time] Tree building completed in {time.time() - start_total:.2f} seconds")

    # save outputs
    #save_tree_output(output_prefix, taxon_nodes, linear_genomes, ordered_nodes)

    # report anomalies
    a = getattr(root_node, "anomalies_total", None)
    if a and (a["n_orphans"] or a["n_self_parents"] or a["n_unreachable"]):
        print("\n=== Anomalies detected in this tree ===")

        def show_nodes(label, tids):
            if not tids:
                return
            print(f" {label}: {len(tids)} node(s)")
            for tid in tids[:5]:  # limit output
                n = taxon_nodes.get(tid)
                if not n:
                    print(f"   • tax_id={tid} (missing node record)")
                    continue
                start, end = getattr(n, "genome_span", (None, None))
                if start is not None and end is not None:
                    count = end - start
                    accs = []
                    if 'linear_genomes' in locals() and linear_genomes:
                        accs = [getattr(g, 'accession', str(g)) for g in linear_genomes[start:start + 3]]
                    acc_str = ", ".join(accs) if accs else "-"
                    print(f"   • {n.name} (tax_id={tid}) — genomes={count}, e.g. {acc_str}")
                else:
                    print(f"   • {n.name} (tax_id={tid}) — no span / no genomes")

        show_nodes("Orphans", a.get("orphans", []))
        show_nodes("Self-parents", a.get("self_parents", []))
        show_nodes("Unreachable", a.get("unreachable", []))
        print("=== End anomalies ===\n")

    return root_node, taxon_nodes, linear_genomes, ordered_nodes, ordered_leaves, name_to_taxids


'''---build tree with phantom nodes---'''
def build_tree_with_phantoms(taxonomy_file, genome_file, output_prefix="output_prefix_phantom", ensembl_tsv_files=None):
    start_total = time.time()

    # Step 1: Parse taxonomy
    t0 = time.time()
    root_node, taxon_nodes = parse_taxonomy(taxonomy_file)
    print(f"[Step 1] Parsed taxonomy in {time.time() - t0:.2f} seconds")

    # Step 2: Attach genomes
    t0 = time.time()
    attach_genomes_to_nodes(genome_file, taxon_nodes)

    # Attach Ensembl genomes
    if ensembl_tsv_files:
        attach_ensembl_tsv_to_nodes(ensembl_tsv_files, taxon_nodes)

    print(f"[Step 2] Attached genomes in {time.time() - t0:.2f} seconds")

    # Step 3: Collect genomes and genome positions
    t0 = time.time()
    linear_genomes, genome_positions = collect_genomes(root_node)
    print(f"[Step 3] Collected genome positions in {time.time() - t0:.2f} seconds")

    # Step 4: Interpolate ranks and assign genome spans
    t0 = time.time()
    annotate_tree_with_ranks_and_spans(root_node, genome_positions)
    print(f"[Step 4] Annotated tree (interpolation & spans) in {time.time() - t0:.2f} seconds")

    # Step 5: Fill in phantom nodes
    t0 = time.time()
    fill_missing_ranks(root_node, taxon_nodes)
    print(f"[Step 5] Inserted phantom nodes in {time.time() - t0:.2f} seconds")

    # Step 6: Re-annotate after phantom changes
    t0 = time.time()
    annotate_tree_with_ranks_and_spans(root_node, genome_positions)
    print(f"[Step 6] Re-annotated tree after phantom insertion in {time.time() - t0:.2f} seconds")

    # Step 7: Linearize tree
    t0 = time.time()
    ordered_leaves, ordered_nodes = linearize_tree(root_node)
    print(f"[Step 8] Linearized tree in {time.time() - t0:.2f} seconds")

    # Step 8: Link siblings
    t0 = time.time()
    link_global_rank_siblings(ordered_nodes)
    print(f"[Step 9] Linked siblings in {time.time() - t0:.2f} seconds")

    # Step 9: Build name index
    t0 = time.time()
    name_to_taxids = defaultdict(list)
    for node in taxon_nodes.values():
        if getattr(node, "phantom", False):
            continue
        orig = getattr(node, "original_name", None)
        if orig:
            name_to_taxids[orig.lower()].append(node.tax_id)
    print(f"[Step 10] Built name index in {time.time() - t0:.2f} seconds")

    print(f"\n[Total time] Tree building completed in {time.time() - start_total:.2f} seconds")

    # save outputs
    #save_tree_output(output_prefix, taxon_nodes, linear_genomes, ordered_nodes)

    # report anomalies
    a = getattr(root_node, "anomalies_total", None)
    if a and (a["n_orphans"] or a["n_self_parents"] or a["n_unreachable"]):
        print("\n=== Anomalies detected in this tree ===")

        def show_nodes(label, tids):
            if not tids:
                return
            print(f" {label}: {len(tids)} node(s)")
            for tid in tids[:5]:  # limit output
                n = taxon_nodes.get(tid)
                if not n:
                    print(f"   • tax_id={tid} (missing node record)")
                    continue
                start, end = getattr(n, "genome_span", (None, None))
                if start is not None and end is not None:
                    count = end - start
                    accs = []
                    if 'linear_genomes' in locals() and linear_genomes:
                        accs = [getattr(g, 'accession', str(g)) for g in linear_genomes[start:start + 3]]
                    acc_str = ", ".join(accs) if accs else "-"
                    print(f"   • {n.name} (tax_id={tid}) — genomes={count}, e.g. {acc_str}")
                else:
                    print(f"   • {n.name} (tax_id={tid}) — no span / no genomes")

        show_nodes("Orphans", a.get("orphans", []))
        show_nodes("Self-parents", a.get("self_parents", []))
        show_nodes("Unreachable", a.get("unreachable", []))
        print("=== End anomalies ===\n")

    return root_node, taxon_nodes, linear_genomes, ordered_nodes, ordered_leaves, name_to_taxids



if __name__ == "__main__":
    taxonomy_input = "taxonomy_all_new.jsonl"
    genomes_input = "genome_all_new.jsonl"
    output_prefix_basic = "tree_basic"
    output_prefix_phantom = "tree_with_phantoms"
    ensembl_txts = [
        "species_EnsemblBacteria.txt",
        "species_EnsemblFungi.txt",
        "species_EnsemblMetazoa.txt",
        "species_EnsemblPlants.txt",
        "species_EnsemblProtists.txt",
        "species_EnsemblVertebrates.txt",
    ]

    print("=== Building Basic Taxonomic Tree (No Phantom Nodes) ===\n")


    build_tree_basic(
        taxonomy_file=taxonomy_input,
        genome_file=genomes_input,
        output_prefix=output_prefix_basic,
        ensembl_tsv_files=ensembl_txts
    )

    # print("\n=== Basic tree building completed. Output saved. ===")



    print("=== Building Taxonomic Tree with Phantom Nodes ===\n")

    build_tree_with_phantoms(
        taxonomy_file=taxonomy_input,
        genome_file=genomes_input,
        output_prefix=output_prefix_phantom,
        ensembl_tsv_files=ensembl_txts
    )

    #print("\n=== Tree with phantom nodes building completed. Output saved. ===")
