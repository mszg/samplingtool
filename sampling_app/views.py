from django.shortcuts import render
from django.http import JsonResponse, HttpResponse
from django.conf import settings
import os
import datetime
import random
import string
import re

# --------------------------
# Helper functions
# --------------------------

def generate_job_id():
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S-")
    random_string = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
    return timestamp + random_string

def clean_filename(filename):
    filename = filename.strip()
    return re.sub(r'[ ,;:<>\$%/\|\\()\[\]{}&*?!\'"\`]', '_', filename)

def create_email_directory(email="guest"):
    safe_email = email.replace("@", "_at_").replace(".", "_dot_")
    email_dir = os.path.join(settings.MEDIA_ROOT, safe_email)
    os.makedirs(email_dir, exist_ok=True)
    return email_dir

# --------------------------
# Views
# --------------------------

def home(request):
    return render(request, "sampling_app/home.html")

def run_sampling(request):
    if request.method != "POST":
        return HttpResponse("Invalid request", status=405)

    taxon = request.POST.get("taxon", "").strip()
    rank = request.POST.get("rank", "").strip()
    genomes = request.POST.get("genomes", "").strip()

    if not taxon:
        return JsonResponse({"error": "Taxon is required"}, status=400)

    # Create directory for guest user
    email_dir = create_email_directory("guest")

    job_id = generate_job_id()
    job_dir = os.path.join(email_dir, job_id)
    os.makedirs(job_dir, exist_ok=True)

    # Output file where Helene's tool will write results
    output_tsv = os.path.join(job_dir, "sampling_results.tsv")

    import subprocess
    try:
        # Call Helene's sampling.py
        cmd = [
            "python",
            str(settings.BASE_DIR_BACKEND / "sampling.py"),
            "--tree", "phantom",
            "sample",
            "--query", taxon,
            "--rank", rank,
            "--per_taxon", genomes,
            "--export_tsv", output_tsv,
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            return JsonResponse({
                "error": "Sampling tool failed.",
        	"stderr": result.stderr,
        	"stdout": result.stdout,
        	"cmd": result.args,
            }, status=500)

    except Exception as exc:
        return JsonResponse({"error": str(exc)}, status=500)

    # Build link for browser download
    rel_path = os.path.relpath(output_tsv, settings.MEDIA_ROOT).replace("\\", "/")
    download_url = settings.MEDIA_URL + rel_path

    return JsonResponse({
        "message": "Sampling finished successfully.",
        "download_url": download_url,
        "job_id": job_id,
    })

