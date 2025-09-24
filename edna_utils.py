import math
import os
from Bio.Blast import NCBIWWW, NCBIXML   # Needs Biopython
from Bio import Entrez

# ---------------- API CONFIG ----------------
Entrez.email = "atharvanimbalkarsr4@gmail.com"   # Change to your email
if "NCBI_API_KEY" in os.environ:
    Entrez.api_key = os.environ["d13f178d009574dccc6a3737da9f84c96a08"]

# ---------------- FASTA Parsing ----------------
def parse_fasta(fasta_str):
    sequences = {}
    seq_id = None
    seq = []
    for line in fasta_str.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq_id:
                sequences[seq_id] = "".join(seq)
            seq_id = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
    if seq_id:
        sequences[seq_id] = "".join(seq)
    return sequences

def clean_species_name(name):
    tokens = name.split()
    if len(tokens) >= 2:
        return " ".join(tokens[:2])
    return name

def mock_blast(sequence):
    if "ATG" in sequence:
        return {
            "species": "Bathyporeia spp.",
            "hit_def": "Bathyporeia spp. [mock hit]",
            "accession": "MOCK123"
        }
    return {
        "species": "Unknown sp.",
        "hit_def": "No reliable hit (mock)",
        "accession": None
    }

def run_blast(sequence):
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence, hitlist_size=1)
        blast_record = NCBIXML.read(result_handle)
        if blast_record.alignments:
            first_hit = blast_record.alignments[0]
            hit_def = first_hit.hit_def
            species_name = clean_species_name(hit_def)
            return {
                "species": species_name,
                "hit_def": hit_def,
                "accession": first_hit.accession
            }
    except Exception as e:
        return {
            "species": "Error",
            "hit_def": str(e),
            "accession": None
        }
    return {
        "species": "Unknown sp.",
        "hit_def": "No hit",
        "accession": None
    }

def get_taxonomy_entrez(scientific_name):
    try:
        handle = Entrez.esearch(db="taxonomy", term=scientific_name)
        record = Entrez.read(handle)
        if record["IdList"]:
            taxid = record["IdList"][0]
            handle = Entrez.efetch(db="taxonomy", id=taxid)
            records = Entrez.read(handle)
            lineage = {d['Rank']: d['ScientificName'] for d in records[0]['LineageEx']}
            return {
                "kingdom": lineage.get("kingdom", "Unknown"),
                "phylum": lineage.get("phylum", "Unknown"),
                "class": lineage.get("class", "Unknown")
            }
    except Exception:
        pass
    return {"kingdom": "Unknown", "phylum": "Unknown", "class": "Unknown"}

def shannon_index(counts):
    total = sum(counts.values())
    return -sum((c / total) * math.log(c / total) for c in counts.values() if c > 0)

def simpson_index(counts):
    total = sum(counts.values())
    return 1 - sum((c / total) ** 2 for c in counts.values() if c > 0)

def chao1(counts):
    counts_list = list(counts.values())
    S_obs = sum(1 for c in counts_list if c > 0)
    F1 = sum(1 for c in counts_list if c == 1)
    F2 = sum(1 for c in counts_list if c == 2)
    if F2 == 0:
        return S_obs + (F1 * (F1 - 1)) / 2
    return S_obs + (F1 ** 2) / (2 * F2)
