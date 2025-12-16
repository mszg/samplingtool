from collections import Counter
import random
import time
import bisect
import argparse
from classes_ranks_definition import *
from tree_building import build_tree_with_phantoms, build_tree_basic
import os, csv, subprocess, glob

class HelpFmt(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter): # formatter class for displaying default values in the help output
    pass

'''---functions for input---'''
def find_taxon(query, taxon_nodes, name_to_taxids):
    if isinstance(query, int): # if query is a tax_id (int) it looks up in taxon_nodes
        return [taxon_nodes.get(query)] if query in taxon_nodes else []
    elif isinstance(query, str): # if query is a name (str) it goes to find_taxon_by_name
        return find_taxon_by_name(query, name_to_taxids, taxon_nodes)
    return None

def find_taxon_by_name(name, name_to_taxids, taxon_nodes): # converts name to lowercase and looks tax_id up (in dict created from tree construction), returns TaxonNode
    matches = []
    name_key = name.lower()
    taxid_list = name_to_taxids.get(name_key, [])
    for tid in taxid_list:
        node = taxon_nodes.get(tid)
        if node:
            matches.append(node)
    return matches

def collect_taxids_from_args(args, name_to_taxids): # user input parser: input can be name/s, taxid/s or a file
    names = list(args.name or [])
    taxids = list(args.taxid or [])

    if getattr(args, "file", None): # if file is passed, it is opened and read line by line
        for raw in args.file:
            s = raw.strip() # remove whitespace
            if not s or s.startswith("#"): # empty lines/comments are skipped
                continue
            if s.isdigit(): # a digit is treated as tax_id
                taxids.append(int(s))
            else:
                names.append(s) # else: name

    if not names and not taxids: # if both lists are emtpy now: Error
        print("[ERROR] Provide at least one --name, --taxid, or --file")
        exit(2)

    # lookup, which is a lowercase key version of name_to_taxids dic
    resolved_taxids = []
    for nm in names:
        tid_list = name_to_taxids.get(nm.lower(), []) # check if name has exact match
        if not tid_list:
            msg = f"[WARN] Taxon name not found: '{nm}'"
            if getattr(args, "strict", False):
                print(msg)
                exit(1)
            else:
                print(msg)
                continue
        resolved_taxids.extend(tid_list)  # flatten

    # combine, preserve order and de-duplicate
    all_taxids = resolved_taxids + taxids
    seen, unique = set(), []
    for tid in all_taxids:
        if tid not in seen:
            unique.append(tid)
            seen.add(tid)
    return unique

def get_genomes_for_node(node, linear_genomes): # reads genome_span (tuple, start and end indices from linear_genomes list) from node
    start, end = getattr(node, "genome_span", (None, None))
    if start is not None and end is not None:
        return linear_genomes[start:end] # slices list of all genomes to return only those under node’s subtree
    return []

def collect_exclusions(args, name_to_taxids, taxon_nodes, linear_genomes): # takes input and creates tuples of what to exclude
    names = list(args.exclude_name or [])
    taxids = list(args.exclude_taxid or [])

    if getattr(args, "exclude_file", None):
        for raw in args.exclude_file:
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            if s.isdigit():
                taxids.append(int(s))
            else:
                names.append(s)

    # resolve names to taxids
    for nm in names:
        tid_list = name_to_taxids.get(nm.lower(), [])
        if not tid_list:
            print(f"[WARN] Exclusion not found: '{nm}'")
            continue
        taxids.extend(tid_list)

    # combine, preserve order and de-duplicate
    excluded_taxids = []
    seen = set()
    for tid in taxids:
        if tid not in seen:
            excluded_taxids.append(tid)
            seen.add(tid)

    # build intervals + accessions which should be excluded
    excluded_leaf_intervals = []
    excluded_accessions = set()
    for tid in excluded_taxids:
        node = taxon_nodes.get(tid)
        if not node:
            print(f"[WARN] Exclusion taxid not in tree: {tid}")
            continue

        # leaf span interval, later used to drop subtrees which fall inside interval
        if node.span and node.span[0] is not None and node.span[1] is not None:
            excluded_leaf_intervals.append((node.span[0], node.span[1]))

        # genomes under this subtree: excluded accession set
        for g in get_genomes_for_node(node, linear_genomes):
            excluded_accessions.add(g.accession)

    return excluded_taxids, excluded_leaf_intervals, excluded_accessions # list of taxid, list of tuples (start_leaf_index, end_leaf_index), set of genome accession strings

def normalize_source_mode(source_mode): # normalize various representations to the three canonical values:
    if not source_mode:
        return "NCBI"
    if source_mode in ("NCBI", "Ensembl", "both"):
        return source_mode
    s = str(source_mode).strip().lower()
    if s == "ncbi":
        return "NCBI"
    if s == "ensembl":
        return "Ensembl"
    if s == "both":
        return "both"
    # fallback
    return "NCBI"

'''---functions for information prints---'''
def print_tree(node, depth=0, max_depth=3):
    indent = "  " * depth
    print(f"{indent}{node.name} (tax_id={node.tax_id}, rank={node.rank})")

    if max_depth == 0 or depth < max_depth: # if max_depth is 0 than print whole subtree
        for child in node.children:
            print_tree(child, depth + 1, max_depth)

def print_taxon_node_info(node):
    print(f"\n--- TaxonNode Info ---")
    print(f"Name: {node.name}")
    print(f"Tax ID: {node.tax_id}")
    print(f"Parent ID: {node.parent_id}")
    print(f"Rank: {node.rank}")
    print(f"Original Rank: {node.original_rank}")
    print(f"Rank Level: {node.rank_level}")
    print(f"Phantom: {node.phantom}")
    print(f"Interpolated: {node.interpolated}")
    print(f"Span (leaves): {node.span}")
    print(f"Sibling: {node.sibling.name}" )
    print(f"Genome Span (genomes): {getattr(node, 'genome_span', (None, None))}")
    print(f"Nearest Descendant Rank: {getattr(node, 'nearest_descendant_rank', None)}")
    print(f"# Children: {len(node.children)}")
    print(f"# Genomes: {len(node.genomes)}")

    if node.genomes:
        print("\nExample Genomes:")
        for g in node.genomes[:3]:  # show up to 3
            ref_flag = "Yes" if g.is_reference else "No"
            print(f" - {g.name} (accession={g.accession}, refseq={ref_flag})")

def print_lineage_to_root(taxon_id, taxon_nodes):
    path = []
    current = taxon_nodes.get(taxon_id)

    while current:
        path.append(current)
        if current.tax_id == current.parent_id:
            break  # reached root
        current = taxon_nodes.get(current.parent_id)

    print("\n--- Lineage to Root (node to root) ---")
    for node in path:
        print(f"Tax ID: {node.tax_id}, Name: {node.name}, Rank: {node.rank.name}")

def print_siblings(node):
    print(f"\nSiblings for node: {node.name} (rank={node.rank})")
    sibling = node.sibling
    while sibling:
        print(f" -> {sibling.name} (tax_id={sibling.tax_id}, rank={sibling.rank})")
        sibling = sibling.sibling


'''---functions for searching of nodes---'''
def select_nodes(node, target_rank, method, ordered_nodes, excluded_leaf_intervals=None): # return list of (node, is_fallback) like current sampler’s discovery phase
    if target_rank in RANK_ORDER:
        rank_attr = "rank"
        target_rank_level = RANK_ORDER.index(target_rank)
    else:
        rank_attr = "original_rank"
        target_rank_level = None

    matched_nodes = []

    if method == "DFS":  # closest descendants of query node that match the target_rank, using DFS traversal, allow fallback to higher ranks if needed
        matched_nodes_raw = get_closest_descendants_with_genomes(node, target_rank_level, node.span)  # finds child nodes at the desired rank or fallback if none
        matched_nodes = [(n, is_fallback) for n, is_fallback in matched_nodes_raw]  # returns list of (node, is_fallback) pairs

    elif method == "list":  # iteration over all nodes (ordered_nodes) to find those inside the query's span with matching rank
        start_end = getattr(node, "span", (None, None))
        if start_end is None or start_end[0] is None or start_end[1] is None:
            return [], target_rank_level, rank_attr
        start, end = start_end

        for n in ordered_nodes:
            if (getattr(n, rank_attr) == target_rank and n.span and n.span[0] is not None and n.span[
                1] is not None and start <= n.span[0] < n.span[1] <= end and n.genome_span != (None, None)):
                matched_nodes.append((n, False))

    elif method == "bisect":  # span lookup using binary search over sorted nodes
        start_end = getattr(node, "span", (None, None))
        if start_end is None or start_end[0] is None or start_end[1] is None:
            return [], target_rank_level, rank_attr
        start, end = start_end

        def get_span_start(node):
            return node.span[0] if node.span else -1

        i = bisect.bisect_left(ordered_nodes, start, key=get_span_start)  # find first index in ordered_nodes where span[0] is bigger/equal to start
        j = bisect.bisect_right(ordered_nodes, end - 1, key=get_span_start)  # like left but for upper bound
        span_nodes = ordered_nodes[i:j]

        matched_nodes = [
            (n, False)
            for n in span_nodes
            if (getattr(n, rank_attr) == target_rank and n.span and n.span[0] is not None and n.span[
                1] is not None and start <= n.span[0] < n.span[1] <= end and  n.genome_span != (None, None))
        ]

    elif method == "sibling":  # sample nodes of given rank linked via .sibling attribute
        start_end = getattr(node, "span", (None, None))
        if start_end is None or start_end[0] is None or start_end[1] is None:
            return [], target_rank_level, rank_attr
        start, end = start_end

        first = None
        for n in ordered_nodes:
            if getattr(n, rank_attr) == target_rank and n.span and start <= n.span[0] < n.span[1] <= end:
                first = n  # start from first match inside query taxon subtree
                break

        if first:
            current = first # follow sibling chain, accepting all same-rank siblings even if gaps exist
            while (current
                    and current.span
                    and start <= current.span[0] < current.span[1] <= end # validate that current is still part of original subtree
                    and getattr(current, rank_attr) == target_rank):
                    if current.genome_span != (None, None):
                        matched_nodes.append((current, False))
                    current = current.sibling

    else:
        # Unknown method -> no matches
        return [], target_rank_level, rank_attr

    if excluded_leaf_intervals:
        def is_inside_excluded(n):
            if not (n.span and n.span[0] is not None and n.span[1] is not None):
                return False
            s, e = n.span
            for xs, xe in excluded_leaf_intervals:
                if xs <= s and e <= xe:
                    return True
            return False

        matched_nodes = [(n, fb) for (n, fb) in matched_nodes if not is_inside_excluded(n)]

    return matched_nodes, target_rank_level, rank_attr


'''---functions for searching and sampling---'''
def get_closest_descendants_with_genomes(node, target_rank_level, query_span): # finds closest descendant nodes that matches target rank level
    result = [] # stores tuples (node, is_fallback)
    seen_fb = set()

    def dfs(n): # postorder
        if not n.span or not (query_span[0] <= n.span[0] < n.span[1] <= query_span[1]): # ensures the node lies within the specified span
            return ('none', [])

        has_genomes = (n.genome_span != (None, None))

        if n.rank_level == target_rank_level and has_genomes: # if node exactly at desired rank -> added to result and marked as exact
            result.append((n, False))
            return ('exact', [])

        # explore children
        child_cands = []  # accumulate per-child closest fallbacks
        saw_exact = False

        for child in n.children:
            status, cands = dfs(child)
            if status == 'exact': # if exact match is found below -> prefer that
                saw_exact = True
            elif status == 'cands':# cands is already the closest level for that child’s branch
                child_cands.extend(cands)

        if saw_exact: # use exact over fallback
            for c in child_cands:
                if c not in seen_fb:
                    result.append((c, True))
                    seen_fb.add(c)
            return ('exact', [])

        if (n.rank_level is not None and target_rank_level is not None
                and n.rank_level > target_rank_level and has_genomes):
            return ('cands', [n])

        if child_cands:
            return ('cands', child_cands)

        return ('none', [])

    # if it returns candidates (no exact anywhere under some branches) it records all nodes at the closest fallback level
    status, cands = dfs(node)
    if status == 'cands':
        for c in cands:
            if c not in seen_fb:
                result.append((c, True))
                seen_fb.add(c)

    return result

def select_genomes(matched_nodes, per_taxon, linear_genomes,
                   source_mode="NCBI",
                   prefer_reference=False, prefer_higher_level=False,
                   min_assembly_level=None, excluded_accessions=None):

    source_mode = normalize_source_mode(source_mode)

    assembly_priority = {
        "COMPLETE GENOME": 4,
        "CHROMOSOME": 3,
        "SCAFFOLD": 2,
        "CONTIG": 1,
    }

    def level_score(g):
        lvl = (g.assembly_level or "").upper()
        return assembly_priority.get(lvl, 0)

    selected = []
    seen_accessions = set()
    node_to_genomes = {}
    rank_usage = Counter()
    fallback_usage = Counter()
    total_available_genomes = 0

    for n, is_fallback in matched_nodes:
        if source_mode == "Ensembl":
            # gather all Ensembl genomes in the subtree of n
            stack = [n]
            all_genomes = []
            while stack:
                cur = stack.pop()
                g_list = getattr(cur, "genomes", []) or []
                for g in g_list:
                    if getattr(g, "origin_source", None) == "Ensembl":
                        all_genomes.append(g)
                stack.extend(cur.children)
        else:
            all_genomes = get_genomes_for_node(n, linear_genomes)

        total_available_genomes += len(all_genomes)

        # remove accessions already used elsewhere
        filtered = [g for g in all_genomes if g.accession not in seen_accessions]

        # remove globally excluded accessions
        if excluded_accessions:
            filtered = [g for g in filtered if g.accession not in excluded_accessions]

        # Ensembl mode: ignore min_assembly_level and preference flags
        if source_mode != "Ensembl" and min_assembly_level:
            min_level_score = assembly_priority.get(min_assembly_level.upper(), 0)
            filtered = [g for g in filtered if level_score(g) >= min_level_score]

        if not filtered:
            continue

        count = min(per_taxon, len(filtered))
        if count <= 0:
            continue

        # Ensembl mode: simple random sampling per node
        if source_mode == "Ensembl":
            chosen = random.sample(filtered, count)

        # prefer_reference (NCBI modes only)
        elif prefer_reference:
            refs = [g for g in filtered if getattr(g, "is_reference", False)]
            non_refs = [g for g in filtered if not getattr(g, "is_reference", False)]

            chosen = []
            remaining = count

            # 1) use reference genomes first
            if refs:
                take = min(remaining, len(refs))
                chosen.extend(random.sample(refs, take))
                remaining -= take

            # 2) if still need more, fill from best assembly levels
            if remaining > 0 and non_refs:
                levels = sorted({level_score(g) for g in non_refs}, reverse=True)
                used = {g.accession for g in chosen}
                for lvl in levels:
                    if remaining <= 0:
                        break
                    pool_lvl = [g for g in non_refs if level_score(g) == lvl and g.accession not in used]
                    if not pool_lvl:
                        continue
                    take = min(remaining, len(pool_lvl))
                    picked = random.sample(pool_lvl, take)
                    chosen.extend(picked)
                    used.update(g.accession for g in picked)
                    remaining -= take

        # prefer_higher_level only (NCBI modes only)
        elif prefer_higher_level:
            chosen = []
            remaining = count
            levels = sorted({level_score(g) for g in filtered}, reverse=True)
            used = set()

            for lvl in levels:
                if remaining <= 0:
                    break
                pool_lvl = [g for g in filtered if level_score(g) == lvl and g.accession not in used]
                if not pool_lvl:
                    continue
                take = min(remaining, len(pool_lvl))
                picked = random.sample(pool_lvl, take)
                chosen.extend(picked)
                used.update(g.accession for g in picked)
                remaining -= take

        # NCBI, no preferences: random
        else:
            chosen = random.sample(filtered, count)

        # if picked fewer than 'count' but enough genomes exist, top up randomly from remaining filtered
        if len(chosen) < count and len(filtered) >= count:
            remaining_needed = count - len(chosen)
            already = {g.accession for g in chosen}
            leftover = [g for g in filtered if g.accession not in already]
            if leftover:
                extra = random.sample(leftover, min(remaining_needed, len(leftover)))
                chosen.extend(extra)

        # bookkeeping
        for g in chosen:
            seen_accessions.add(g.accession)
        selected.extend(chosen)
        node_to_genomes[n.name] = (n, chosen)

        rank_usage[n.rank.name] += 1
        if is_fallback:
            fallback_usage[n.rank.name] += 1

    return selected, node_to_genomes, rank_usage, fallback_usage, total_available_genomes


'''---functions for preview and print---'''
def count_ensembl_genomes_in_subtree(node):
    stack = [node]
    count = 0
    while stack:
        cur = stack.pop()
        for g in getattr(cur, "genomes", []) or []:
            if getattr(g, "origin_source", None) == "Ensembl":
                count += 1
        stack.extend(cur.children)
    return count

def compute_subtree_counts(node): # returns span and count info for query node
    # leaf span: leaves
    leaves_start, leaves_end = node.span if node.span else (None, None)
    leaves_count = (max(0, leaves_end - leaves_start)
                    if (leaves_start is not None and leaves_end is not None) else 0)

    # genome span: genomes
    genome_start, genome_end = getattr(node, "genome_span", (None, None))
    genomes_under_query = (max(0, genome_end - genome_start)
                           if (genome_start is not None and genome_end is not None) else 0)

    return (leaves_start, leaves_end, leaves_count, genome_start, genome_end, genomes_under_query)

def build_sampling_plan(query, node, target_rank, per_taxon, method, *,
                        source_mode="NCBI",
                        prefer_reference=False, prefer_higher_level=False,
                        min_assembly_level=None, ordered_nodes=None, linear_genomes=None,
                        excluded_leaf_intervals=None, excluded_accessions=None):
    t0 = time.time()
    source_mode = normalize_source_mode(source_mode)

    # normalize excluded_accessions to a set
    effective_excluded_accessions = set(excluded_accessions or [])

    # handle 'both' mode: restrict to NCBI genomes whose accessions are also in Ensembl
    if source_mode == "both":
        # root_node is global, set when building the tree
        stats = getattr(root_node, "genome_source_stats", None)
        if stats is None:
            print("[WARN] source_mode 'both' requested but genome_source_stats "
                "not found on root; falling back to NCBI-only.")
        else:
            ncbi_accs = stats.get("ncbi_accessions", set())
            shared_accs = stats.get("shared_accessions", set())
            extra_excluded = ncbi_accs - shared_accs
            effective_excluded_accessions |= extra_excluded

    # nodes searching with subtree exclusion
    matched_nodes, target_rank_level, rank_attr = select_nodes(node, target_rank, method, ordered_nodes, excluded_leaf_intervals=excluded_leaf_intervals)

    # counts derived directly from discovery result (method-dependent)
    candidate_exact_nodes = sum(1 for _, fb in matched_nodes if not fb)
    candidate_fallback_nodes = sum(1 for _, fb in matched_nodes if fb)

    # spans & genome count (method-independent)
    (leaves_start, leaves_end, leaves_count,
     genome_start, genome_end, genomes_under_query) = compute_subtree_counts(node)

    # selection
    (selected, node_to_genomes, rank_usage, fallback_usage,
     total_available_genomes_in_matched) = select_genomes(
        matched_nodes, per_taxon, linear_genomes,
        source_mode=source_mode,
        prefer_reference=prefer_reference,
        prefer_higher_level=prefer_higher_level,
        min_assembly_level=min_assembly_level,
        excluded_accessions=effective_excluded_accessions)

    plan = {
        "query": query,
        "node": node,
        "target_rank": target_rank,
        "method": method,
        "source_mode": source_mode,
        "matched_nodes": matched_nodes, # list [(node, is_fallback)]
        "target_rank_level": target_rank_level,

        # method-dependent counts straight from select_nodes
        "candidate_exact_nodes": candidate_exact_nodes,
        "candidate_fallback_nodes": candidate_fallback_nodes,

        "excluded_nodes": len(excluded_leaf_intervals) if excluded_leaf_intervals else 0,
        "excluded_accessions": len(excluded_accessions) if excluded_accessions else 0,

        # spans & genome count
        "leaves_start": leaves_start,
        "leaves_end": leaves_end,
        "leaves_count": leaves_count,
        "genome_start": genome_start,
        "genome_end": genome_end,
        "genomes_under_query": genomes_under_query,

        "total_available_genomes_in_matched": total_available_genomes_in_matched,
        "selected": selected, # list[Genomes]
        "node_to_genomes": node_to_genomes, # dict[name]: (node, [Genome])
        "rank_usage": rank_usage, # Counter by rank name
        "fallback_usage": fallback_usage, # Counter for fallback nodes

        "elapsed_sec": time.time() - t0}
    return plan

def print_preview_report(plan): # prints compact summery shown before results
    print("\n=== Sampling Preview ===")
    print(
        f"Query: '{plan['query']}'  | Target rank: {plan['target_rank'].name}  " f"| Method: {plan['method']}  | Source: {plan.get('source_mode', 'NCBI')}")

    if plan["leaves_start"] is None:
        print("- Leaves span: unavailable for this node (no span).")
    else:
        print(f"- Leaves span:  start={plan['leaves_start']}, end={plan['leaves_end']}, leaves={plan['leaves_count']}")

    source = plan.get("source_mode", "NCBI")

    if source == "Ensembl":
        # Ensembl genome count instead of NCBI genome span
        node = plan["node"]  # sampling point
        ens_count = count_ensembl_genomes_in_subtree(node)

        print(f"- Ensembl genomes in subtree: {ens_count}")
    else:
        # NCBI / both modes: keep the normal NCBI-based spans
        if plan["genome_start"] is None:
            print("- Genome span: unavailable for this node (no genome span).")
        else:
            print(
                f"- Genome span:  start={plan['genome_start']}, end={plan['genome_end']}, genomes={plan['genomes_under_query']}")

    print(f"- Candidate nodes to sample from: {len(plan['matched_nodes'])} "
          f"(exact={plan['candidate_exact_nodes']}, fallback={plan['candidate_fallback_nodes']})")

    if plan.get("excluded_nodes", 0) or plan.get("excluded_accessions", 0):
        print(f"- Exclusions: {plan.get('excluded_nodes', 0)} subtree(s), "
              f"{plan.get('excluded_accessions', 0)} genome(s) removed pre-selection")

    print(f"- Genomes that would be selected now: {len(plan['selected'])}")
    print("───────────────────────────────────────────────────────")
    print("Proceed?  [y]es to print results  |  [n]o to change inputs  |  [s]top to exit")


def build_and_preview(q_query, q_rank, q_per_taxon, q_method, q_source_mode, q_min_level, q_pref_ref, q_pref_high, excluded_leaf_intervals=None, excluded_accessions=None):
    # allow query to be either a name or a numeric taxid
    query_for_lookup = q_query
    if isinstance(q_query, str) and q_query.isdigit():
        query_for_lookup = int(q_query)
    nodes = find_taxon(query_for_lookup, taxon_nodes, name_to_taxids)

    if not nodes:
        print(f"[ERROR] Taxon '{q_query}' not found.")
        return None, None
    if len(nodes) > 1:
        print(f"[INFO] Multiple taxa found for name '{q_query}':")
        for i, node in enumerate(nodes):
            print(f"  [{i}] {node.name} (tax_id={node.tax_id}, rank={node.rank.name})")
        try:
            choice = int(input("Select index of taxon to use: ").strip())
            node = nodes[choice]
        except (ValueError, IndexError):
            print("[ERROR] Invalid selection.")
            return None, None
    else:
        node = nodes[0]

    plan = build_sampling_plan(query=q_query, node=node, target_rank=q_rank, per_taxon=q_per_taxon, method=q_method, source_mode=q_source_mode,
        prefer_reference=q_pref_ref, prefer_higher_level=q_pref_high, min_assembly_level=q_min_level,
        ordered_nodes=ordered_nodes, linear_genomes=linear_genomes,
        excluded_leaf_intervals=excluded_leaf_intervals,
        excluded_accessions=excluded_accessions)
    print_preview_report(plan)
    return node, plan


'''---functions for result output and export---'''
out_file_handle = None # global
def log(msg="", end="\n"): # print to console and optionally to a text file
    print(msg, end=end)  # always show in console
    if out_file_handle:
        out_file_handle.write(str(msg) + end)

def emit_sampling_results(plan, report_phantoms=True): # print final selection, without recomputation
    matched_nodes = plan["matched_nodes"]
    target_rank_level = plan["target_rank_level"]
    rank_usage = plan["rank_usage"]
    fallback_usage = plan["fallback_usage"]
    node_to_genomes = plan["node_to_genomes"]
    selected = plan["selected"]

    log(f"\nSelected {len(selected)} genomes from {len(matched_nodes)} nodes in {plan['elapsed_sec']:.2f} seconds.")
    if fallback_usage:
        log(f"Fallbacks used in {sum(fallback_usage.values())} nodes.")

    log("\nActual ranks used:")
    for rank in rank_usage:
        fallback_note = " (fallback)" if rank in fallback_usage else ""
        log(f"  {rank}: {rank_usage[rank]} nodes{fallback_note}")

    total_ref = sum(1 for g in selected if getattr(g, "is_reference", False))
    total_nonref = len(selected) - total_ref
    log(f"\nReference status in selection: {total_ref} reference, {total_nonref} non-reference genomes.")

    level_counts = Counter((g.assembly_level or "Unknown") for g in selected)
    log("Assembly levels in selection:")
    for level, cnt in level_counts.most_common():
        log(f"  {level}: {cnt} genomes")

    for node_name, (node, genomes) in node_to_genomes.items():
        label = " [Phantom]" if report_phantoms and getattr(node, "phantom", False) else ""
        interpolation_tag = " (interpolated)" if report_phantoms and getattr(node, "interpolated", False) else ""
        fallback_flag = " [FALLBACK]" if node.rank_level != target_rank_level else ""
        log(f"\n{node.name} ({node.rank.name}){label}{interpolation_tag}{fallback_flag}:")
        for g in genomes:
            ref_flag = "Yes" if g.is_reference else "No"
            log(f"  - {g.name} (accession={g.accession}, refseq={ref_flag}, level={g.assembly_level})")


def export_sampled_genomes(plan, output_tsv, base_dir="/data/genomes", genome_ext=".fna", annot_ext=".gff3"): # make a TSV with columns: Name, Genome_Path, Annotation_Path
   # initially fills *intended* paths (base_dir/<acc>/<acc>.<ext>), correct later to real files after download with datasets
    if not output_tsv:
        return
    selected = plan.get("selected", [])
    os.makedirs(os.path.dirname(output_tsv) or ".", exist_ok=True)
    with open(output_tsv, "w", encoding="utf-8") as out:
        out.write("Name\tGenome_Path\tAnnotation_Path\n")
        for g in selected:
            acc = getattr(g, "accession", None)
            if not acc:
                continue
            genome_path = f"{base_dir}/{acc}/{acc}{genome_ext}"
            annot_path  = f"{base_dir}/{acc}/{acc}{annot_ext}"
            out.write(f"{acc}\t{genome_path}\t{annot_path}\n")
    print(f"[Export] Wrote {len(selected)} entries to {output_tsv}")

def fetch_with_datasets(tsv_path, base_dir="/data/genomes", include="gff3,rna,cds,protein,genome,seq-report", skip_existing=True):
    # download each accession (Name column) with the NCBI Datasets CLI: datasets download genome accession <acc> --include <include> --filename <acc>.zip
    # unzips into base_dir/<acc>/. Skips if a .fna already exists (unless skip_existing=False)
    if not os.path.exists(tsv_path):
        print(f"[ERROR] TSV not found: {tsv_path}")
        return

    with open(tsv_path, newline='', encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc = row["Name"]
            target_dir = os.path.join(base_dir, acc)
            os.makedirs(target_dir, exist_ok=True)

            # skip if have at least one .fna
            if skip_existing and glob.glob(os.path.join(target_dir, "**", "*.fna"), recursive=True):
                print(f"[Skip] {acc}: already present in {target_dir}")
                continue

            print(f"[Fetch] {acc} -> {target_dir}")
            zip_path = os.path.join(target_dir, f"{acc}.zip")
            datasets_exe = "./datasets"
            cmd = [datasets_exe, "download", "genome", "accession", acc,
                "--include", include,
                "--filename", zip_path,
                "--no-progressbar"]
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError:
                print(f"[WARN] Failed to download {acc}")
                continue

            # unzip and cleanup
            subprocess.run(["unzip", "-o", zip_path, "-d", target_dir], check=True)
            os.remove(zip_path)

    print("[Fetch] Completed all downloads.")

def _pick_best_genome_fasta(dir_for_acc: str) -> str | None: # preference order: *_genomic.fna  >  any *.fna
    # datasets layout usually: <acc>/ncbi_dataset/data/<acc>/*_genomic.fna
    candidates = glob.glob(os.path.join(dir_for_acc, "**", "*_genomic.fna"), recursive=True)
    if candidates:
        return sorted(candidates)[0]
    candidates = glob.glob(os.path.join(dir_for_acc, "**", "*.fna"), recursive=True)
    return sorted(candidates)[0] if candidates else None

def _pick_best_annotation(dir_for_acc: str) -> str | None: # prefer *.gff3 if present, otherwise *.gff
    gff3 = glob.glob(os.path.join(dir_for_acc, "**", "*.gff3"), recursive=True)
    if gff3:
        return sorted(gff3)[0]
    gff = glob.glob(os.path.join(dir_for_acc, "**", "*.gff"), recursive=True)
    return sorted(gff)[0] if gff else None

def rewrite_tsv_with_actual_paths(tsv_path, base_dir="/data/genomes"): # locate the real downloaded files under base_dir/<acc>/ncbi_dataset/..., and rewrite the TSV so Genome_Path/Annotation_Path point to real files
    # read rows
    rows = []
    with open(tsv_path, newline='', encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)

    # rewrite paths
    fixed = 0
    for row in rows:
        acc = row["Name"]
        acc_dir = os.path.join(base_dir, acc)
        if not os.path.isdir(acc_dir):
            print(f"[WARN] No directory for {acc} at {acc_dir}")
            continue

        real_fna = _pick_best_genome_fasta(acc_dir)
        real_gff = _pick_best_annotation(acc_dir)

        if real_fna:
            row["Genome_Path"] = real_fna
        else:
            print(f"[WARN] No .fna found for {acc}")

        if real_gff:
            row["Annotation_Path"] = real_gff
        else:
            print(f"[WARN] No .gff/.gff3 found for {acc}")

        if real_fna or real_gff:
            fixed += 1

    # write back
    tmp_out = tsv_path + ".tmp"
    with open(tmp_out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["Name", "Genome_Path", "Annotation_Path"], delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow({
                "Name": row["Name"],
                "Genome_Path": row["Genome_Path"],
                "Annotation_Path": row["Annotation_Path"],
            })
    os.replace(tmp_out, tsv_path)
    print(f"[Export] Rewrote TSV with discovered paths: {tsv_path} (updated {fixed} rows)")



if __name__ == "__main__":
    # option to choose if one wants to sample or get taxon information
    parser = argparse.ArgumentParser(
        description="Work with a taxonomic tree: sample genomes or print node info.",
        formatter_class=HelpFmt,
        epilog=(
            "Examples:\n"
            "  # Preview sampling 2 genomes per FAMILY under Bacteria with DFS in the phantom tree, in an interactive way\n"
            "  python sampling.py --tree phantom sample --query Bacteria --rank FAMILY --per_taxon 2 --method DFS --interactive\n\n"
            "  # Non-interactive sampling with preferences and exclusions\n"
            "  python sampling.py --tree basic sample --query 'Escherichia' --rank SPECIES \\\n"
            "      --per_taxon 3 -- method DFS --prefer_reference --prefer_higher_level \\\n"
            "      --exclude_name 'Escherichia coli' --seed 42\n\n"
            "  # Print node info + lineage and a small subtree\n"
            "  python sampling.py --tree phantom info --name Bacteria --lineage --subtree --subtree_depth 2\n"
        )
    )
    parser.add_argument("--tree", choices=["phantom", "basic"], required=True, help="Which input tree to use: 'phantom' or 'basic'")
    parser.add_argument("--out", type=str, help="Optional path to also write output to a plain text file")

    subparsers = parser.add_subparsers(dest="command", required=True)
    p_sample = subparsers.add_parser("sample", description="Sample genomes by rank from a taxonomic tree.")

    p_sample.add_argument("--query", required=True, help="Name or taxid of taxon to sample under (e.g., 'Primates' or '9443')",)
    p_sample.add_argument("--rank", type=str, required=True, help="Target rank to sample from (e.g., GENUS, FAMILY)")
    p_sample.add_argument("--per_taxon", type=int, default=1, help="Number of genomes per taxon to sample")
    p_sample.add_argument("--source_mode",choices=["NCBI", "Ensembl", "both"], default="NCBI",help=("Which genomes to sample from: "
                                                                                                    "'NCBI' (NCBI only; current behaviour), "
                                                                                                    "'Ensembl' (Ensembl genomes only, simple random per node), "
                                                                                                    "'both' (NCBI genomes whose accessions also exist in Ensembl)."),)
    p_sample.add_argument("--method", choices=["DFS", "list", "bisect", "sibling"], default="DFS",  help="Sampling method to use")
    p_sample.add_argument("--prefer_reference", action="store_true")
    p_sample.add_argument("--prefer_higher_level", action="store_true")
    p_sample.add_argument("--min_assembly_level", type=str, help="e.g., CONTIG")
    p_sample.add_argument("--interactive", action="store_true", help="Preview counts first and interactively confirm or change inputs")
    p_sample.add_argument("--seed", type=int, help="Random seed for reproducible sampling")
    p_sample.add_argument("--exclude_name", nargs="+", help="One or more taxon names to exclude (e.g., --exclude_name 'Homo sapiens' Homo)")
    p_sample.add_argument("--exclude_taxid", nargs="+", type=int, help="One or more taxon IDs to exclude")
    p_sample.add_argument("--exclude_file", type=argparse.FileType("r"), help="Optional file with one name or taxid per line (# comments ok)")

    # export + fetch flags
    p_sample.add_argument("--export_tsv", type=str,
                          help="Write selected genomes to this TSV (Name, Genome_Path, Annotation_Path)")
    p_sample.add_argument("--base_dir", type=str, default="/data/genomes",
                          help="Base dir where <acc>/ will be created and data stored")

    # NCBI Datasets integration
    p_sample.add_argument("--fetch_datasets", action="store_true",
                          help="Download selected genomes with the NCBI Datasets CLI")
    p_sample.add_argument("--datasets_include", type=str,
                          default="gff3,rna,cds,protein,genome,seq-report",
                          help="Comma list for NCBI Datasets --include (default: gff3,rna,cds,protein,genome,seq-report)")

    # info subcommand
    p_info = subparsers.add_parser("info", help="Print detailed information for one or more taxon nodes", description="Print taxon node information.")
    p_info.add_argument("--name", nargs="+", help="One or more taxon names (e.g., --name Primates Canidae)")
    p_info.add_argument("--taxid", nargs="+", type=int, help="One or more numeric taxon IDs (e.g., --taxid 9443 9604)")
    p_info.add_argument("--file", type=argparse.FileType("r"), help="Optional file with one name or taxid per line (lines starting with # are ignored)")
    p_info.add_argument("--strict", action="store_true", help="Fail immediately on the first missing/unknown taxon (default: skip and continue)")
    p_info.add_argument("--lineage", action="store_true", help="Also print lineage to root for each taxon")
    p_info.add_argument("--siblings", action="store_true", help="Also print the forward siblings chain for each taxon")
    p_info.add_argument("--subtree", action="store_true", help="Also print the subtree under each taxon")
    p_info.add_argument("--subtree_depth", type=int, default=3, help="Max depth to print for --subtree (default: 3)")

    # linage subcommand
    p_lineage = subparsers.add_parser("lineage", help="Print lineage to root for one or more taxa", description="Print lineage to root.")
    p_lineage.add_argument("--name", nargs="+", help="One or more taxon names")
    p_lineage.add_argument("--taxid", nargs="+", type=int, help="One or more numeric taxon IDs")
    p_lineage.add_argument("--file", type=argparse.FileType("r"), help="Optional file with one name or taxid per line (# comments ok)")
    p_lineage.add_argument("--strict", action="store_true", help="Fail on first missing/unknown taxon (default: warn and continue)")

    # subtree subcommand
    p_subtree = subparsers.add_parser("subtree", help="Print subtree under one or more taxa", description="Print the tree below the selected node(s).")
    p_subtree.add_argument("--name", nargs="+", help="One or more taxon names")
    p_subtree.add_argument("--taxid", nargs="+", type=int, help="One or more numeric taxon IDs")
    p_subtree.add_argument("--file", type=argparse.FileType("r"), help="Optional file with one name or taxid per line (# comments ok)")
    p_subtree.add_argument("--strict", action="store_true", help="Fail on first missing/unknown taxon (default: warn and continue)")
    p_subtree.add_argument("--depth", type=int, default=3, help="Max depth to print (default: 3)")

    args = parser.parse_args()

    if args.out:
        out_file_handle = open(args.out, "w", encoding="utf-8")

    if getattr(args, "seed", None) is not None:
        random.seed(args.seed)

    if getattr(args, "command", None) == "sample" and getattr(args, "source_mode", "NCBI") == "Ensembl":
        if args.prefer_reference or args.prefer_higher_level or args.min_assembly_level:
            print(
                "[WARN] --source_mode Ensembl: "
                "--prefer_reference, --prefer_higher_level and --min_assembly_level "
                "are ignored (Ensembl genomes are treated as uniformly high-quality)."
            )
    total_start = time.time()

    # determine tree prefix
    print(f"\n[INFO] Using tree: {args.tree}")
    prefix = "tree_with_phantoms" if args.tree == "phantom" else "tree_basic"
    ensembl_txts = ["species_EnsemblBacteria.txt",
        "species_EnsemblFungi.txt",
        "species_EnsemblMetazoa.txt",
        "species_EnsemblPlants.txt",
        "species_EnsemblProtists.txt",
        "species_EnsemblVertebrates.txt"]
    # load data
    print(" Loading data...")
    t0 = time.time()
    print("[BUILD] Building tree in memory...")
    if args.tree == "phantom":
        root_node, taxon_nodes, linear_genomes, ordered_nodes, ordered_leaves, name_to_taxids = build_tree_with_phantoms(
            taxonomy_file="taxonomy_all_new.jsonl",
            genome_file="genome_all_new.jsonl",
            output_prefix="tree_with_phantoms",
            ensembl_tsv_files=ensembl_txts)
    else:
        root_node, taxon_nodes, linear_genomes, ordered_nodes, ordered_leaves, name_to_taxids = build_tree_basic(
            taxonomy_file="taxonomy_all_new.jsonl",
            genome_file="genome_all_new.jsonl",
            output_prefix="tree_basic",
            ensembl_tsv_files=ensembl_txts)
    print(f"[Done] Data loaded in {time.time() - t0:.2f} sec")

    # branch by command
    if args.command == "info": # resolve nodes by name or taxid
        taxid_list = collect_taxids_from_args(args, name_to_taxids)

        successes = 0
        for tid in taxid_list:
            node = taxon_nodes.get(tid)
            if node is None:
                msg = f"[WARN] Taxon ID not found in tree: {tid}"
                if args.strict:
                    print(msg)
                    exit(1)
                else:
                    print(msg)
                    continue

            print_taxon_node_info(node)

            if getattr(args, "siblings", False):
                print_siblings(node)

            if getattr(args, "lineage", False):
                print_lineage_to_root(tid, taxon_nodes)

            if getattr(args, "subtree", False):
                print(f"\n--- Subtree (max_depth={args.subtree_depth}) for {node.name} [{node.tax_id}] ---")
                print_tree(node, depth=0, max_depth=args.subtree_depth)

            successes += 1

        print(f"\n[INFO] Completed info for {successes} of {len(taxid_list)} requested taxa.")
        print(f"\n[Total Time] All steps completed in {time.time() - total_start:.2f} sec")

    elif args.command == "lineage":
        taxid_list = collect_taxids_from_args(args, name_to_taxids)

        successes = 0
        for tid in taxid_list:
            node = taxon_nodes.get(tid)
            if node is None:
                msg = f"[WARN] Taxon ID not found in tree: {tid}"
                if args.strict:
                    print(msg)
                    exit(1)
                else:
                    print(msg)
                    continue

            print_lineage_to_root(tid, taxon_nodes)
            successes += 1

        print(f"\n[INFO] Completed lineage for {successes} of {len(taxid_list)} requested taxa.")
        print(f"\n[Total Time] All steps completed in {time.time() - total_start:.2f} sec")

    elif args.command == "subtree":
        taxid_list = collect_taxids_from_args(args, name_to_taxids)

        successes = 0
        for tid in taxid_list:
            node = taxon_nodes.get(tid)
            if node is None:
                msg = f"[WARN] Taxon ID not found in tree: {tid}"
                if args.strict:
                    print(msg)
                    exit(1)
                else:
                    print(msg)
                    continue

            depth_info = "entire tree" if args.depth == 0 else f"max_depth={args.depth}"
            print(f"\n--- Subtree ({depth_info}) for {node.name} [{node.tax_id}] ---")
            print_tree(node, depth=0, max_depth=args.depth)
            successes += 1

        print(f"\n[INFO] Completed subtree for {successes} of {len(taxid_list)} requested taxa.")
        print(f"\n[Total Time] All steps completed in {time.time() - total_start:.2f} sec")

    elif args.command == "sample":
        # convert rank to enum once up front
        try:
            target_rank = RankType[args.rank.upper()]
        except KeyError:
            print(f"[ERROR] Invalid rank: '{args.rank}'. Must be one of: {[r.name for r in RankType]}")
            exit(1)

        # working copy of parameters that the user can change in-loop
        q_query = args.query
        q_rank = target_rank
        q_per_taxon = args.per_taxon
        q_method = args.method
        q_source_mode = args.source_mode
        q_min_level = args.min_assembly_level
        q_pref_ref = args.prefer_reference
        q_pref_high = args.prefer_higher_level
        report_phantoms = (args.tree == "phantom")
        q_exclude_name = list(args.exclude_name or [])
        q_exclude_taxid = list(args.exclude_taxid or [])
        q_exclude_file_path = args.exclude_file.name if getattr(args, "exclude_file", None) else None
        q_seed = args.seed  # may be None

        excluded_taxids = []
        excluded_leaf_intervals = []
        excluded_accessions = set()
        if any(getattr(args, k, None) for k in ("exclude_name", "exclude_taxid", "exclude_file")):
            excluded_taxids, excluded_leaf_intervals, excluded_accessions = collect_exclusions(
                args, name_to_taxids, taxon_nodes, linear_genomes)
        if not args.interactive:
            # build plan and immediately emit results
            node, plan = build_and_preview(q_query, q_rank, q_per_taxon, q_method, q_source_mode,
                q_min_level, q_pref_ref, q_pref_high,
                excluded_leaf_intervals=excluded_leaf_intervals,
                excluded_accessions=excluded_accessions)
            if plan is None:
                exit(1)
            emit_sampling_results(plan, report_phantoms=report_phantoms)
            # 1) export a TSV (intended paths first)
            if getattr(args, "export_tsv", None):
                export_sampled_genomes(
                    plan,
                    args.export_tsv,
                    base_dir=args.base_dir,
                    genome_ext=None,
                    annot_ext=None)

            print(f"\n[Total Time] All steps completed in {time.time() - total_start:.2f} sec")
            exit(0)

        # interactive loop
        while True:
            if q_seed is not None: # reseed if the user set/changed it
                random.seed(q_seed)

            # rebuild exclusions from current in-loop state (names/taxids/file)
            if (q_exclude_name or q_exclude_taxid or q_exclude_file_path):
                tmp_args = argparse.Namespace(
                    exclude_name=q_exclude_name,
                    exclude_taxid=q_exclude_taxid,
                    exclude_file=open(q_exclude_file_path, "r") if q_exclude_file_path else None)
                _, excluded_leaf_intervals, excluded_accessions = collect_exclusions(tmp_args, name_to_taxids, taxon_nodes, linear_genomes)
                if tmp_args.exclude_file:
                    tmp_args.exclude_file.close()
            else:
                excluded_leaf_intervals = []
                excluded_accessions = set()

            node, plan = build_and_preview(q_query, q_rank, q_per_taxon, q_method, q_source_mode,
                q_min_level, q_pref_ref, q_pref_high,
                excluded_leaf_intervals=excluded_leaf_intervals,
                excluded_accessions=excluded_accessions)
            if plan is None:
                # unresolved query: let user change or stop
                choice = input("[n] change inputs  |  [s] stop: ").strip().lower()
                if choice.startswith("s"):
                    print("\n[INFO] Stopped.")
                    break
                # else continue to edit
            else:
                choice = input("your choice [y/n/s]: ").strip().lower()

            if choice in ("y", "yes"):
                emit_sampling_results(plan, report_phantoms=report_phantoms)
                # 1) export a TSV (intended paths first)
                if getattr(args, "export_tsv", None):
                    export_sampled_genomes(
                        plan,
                        args.export_tsv,
                        base_dir=args.base_dir,
                        genome_ext=".fna",
                        annot_ext=".gff3")

                    # 2) download with NCBI Datasets
                    if getattr(args, "fetch_datasets", False):
                        fetch_with_datasets(
                            args.export_tsv,
                            base_dir=args.base_dir,
                            include=args.datasets_include,
                            skip_existing=True)

                        # 3) rewrite TSV with the real file paths discovered under ncbi_dataset/
                        rewrite_tsv_with_actual_paths(args.export_tsv, base_dir=args.base_dir)
                print(f"\n[Total Time] All steps completed in {time.time() - total_start:.2f} sec")
                break

            if choice in ("s", "stop"):
                print("\n[INFO] Stopped without sampling.")
                print(f"\n[Total Time] All steps completed in {time.time() - total_start:.2f} sec")
                break

            # choice == 'n':  edit parameters interactively
            print("\nEnter new values or press ENTER to keep current in [brackets].")
            new_query = input(f"- query name [{q_query}]: ").strip()
            if new_query:
                q_query = new_query

            new_rank = input(f"- target rank (e.g., GENUS, FAMILY) [{q_rank.name}]: ").strip()
            if new_rank:
                try:
                    q_rank = RankType[new_rank.upper()]
                except KeyError:
                    print(f"  [WARN] Invalid rank '{new_rank}', keeping {q_rank.name}.")

            new_per_taxon = input(f"- genomes per taxon [{q_per_taxon}]: ").strip()
            if new_per_taxon.isdigit():
                q_per_taxon = int(new_per_taxon)

            new_method = input(f"- method (DFS|list|bisect|sibling) [{q_method}]: ").strip()
            if new_method:
                if new_method in ("DFS", "list", "bisect", "sibling"):
                    q_method = new_method
                else:
                    print(f"  [WARN] Unknown method '{new_method}', keeping {q_method}.")

            raw_source = input(f"- source mode (NCBI|Ensembl|both) [{q_source_mode}]: ").strip()
            if raw_source:
                new_mode = normalize_source_mode(raw_source)
                if new_mode != q_source_mode:
                    q_source_mode = new_mode
                    if q_source_mode == "Ensembl" and (q_pref_ref or q_pref_high or q_min_level):
                        print(
                            "[WARN] source_mode=Ensembl: "
                            "prefer_reference, prefer_higher_level and min_assembly_level "
                            "will be ignored; Ensembl genomes are treated as uniformly high quality.")

            new_min = input(f"- min assembly level (CONTIG/… or blank) [{q_min_level}]: ").strip()
            q_min_level = new_min or q_min_level

            new_pref_ref = input(f"- prefer_reference (y/n) [{'y' if q_pref_ref else 'n'}]: ").strip().lower()
            if new_pref_ref in ("y", "yes", "n", "no"):
                q_pref_ref = new_pref_ref.startswith("y")

            new_pref_high = input(f"- prefer_higher_level (y/n) [{'y' if q_pref_high else 'n'}]: ").strip().lower()
            if new_pref_high in ("y", "yes", "n", "no"):
                q_pref_high = new_pref_high.startswith("y")

            print("\n— Exclusions (press ENTER to keep current) —")
            cur_names = ", ".join(q_exclude_name) or "none"
            inp = input(f"- exclude names (comma-separated) [{cur_names}]: ").strip()
            if inp:
                q_exclude_name = [] if inp == "-" else [s.strip() for s in inp.split(",") if s.strip()]

            cur_ids = ", ".join(map(str, q_exclude_taxid)) or "none"
            inp = input(f"- exclude taxids (comma-separated) [{cur_ids}]: ").strip()
            if inp:
                q_exclude_taxid = [] if inp == "-" else [int(s) for s in inp.replace(" ", "").split(",") if
                                                         s.isdigit()]

            cur_file = q_exclude_file_path or "none"
            inp = input(f"- exclude file path [{cur_file}] (use '-' to clear): ").strip()
            if inp:
                q_exclude_file_path = None if inp in ("-", "none") else inp

            print("\n— Random seed (press ENTER to keep current) —")
            cur_seed = "none" if q_seed is None else str(q_seed)
            inp = input(f"- random seed [{cur_seed}] (type 'none' to unset): ").strip().lower()
            if inp:
                if inp in ("none", "-"):
                    q_seed = None
                elif inp.isdigit():
                    q_seed = int(inp)
                else:
                    print(" Seed must be an integer or 'none'; keeping current.")

            # loop back: rebuild plan and show updated preview

    if out_file_handle:
        out_file_handle.close()