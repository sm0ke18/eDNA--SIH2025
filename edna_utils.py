import math
from Bio.Blast import NCBIWWW, NCBIXML   # Needs Biopython
from Bio import Entrez
from ete3 import NCBITaxa



Entrez.email = "atharvanimbalkar18@gmail.com"  # <-- Replace with your email
ncbi = NCBITaxa()

# ---------------- FASTA Parsing ----------------
def parse_fasta(fasta_str):
    """Parse a FASTA string into a dictionary {seq_id: sequence}."""
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
            seq_id = line[1:].split()[0]  # Take first token as ID
            seq = []
        else:
            seq.append(line)
    if seq_id:
        sequences[seq_id] = "".join(seq)
    return sequences

# ---------------- MOCK BLAST ----------------
def mock_blast(sequence):
    """Return mock BLAST results for testing/demo purposes."""
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

# ---------------- REAL BLAST ----------------
def run_blast(sequence):
    """Run BLAST against NCBI nt database and return first hit info."""
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
        blast_record = NCBIXML.read(result_handle)
        if blast_record.alignments:
            first_hit = blast_record.alignments[0]
            return {
                "species": " ".join(first_hit.hit_def.split()[0:2]),
                "hit_def": first_hit.hit_def,
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

# ---------------- Entrez Taxonomy ----------------
def get_taxonomy_entrez(scientific_name):
    """Fetch taxonomy (kingdom, phylum, class) from NCBI Entrez for a given scientific name."""
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

def get_taxonomy_ete(scientific_name):
    """Fetch taxonomy (kingdom, phylum, class) from local ETE NCBI database."""
    try:
        name2taxid = ncbi.get_name_translator([scientific_name])
        if not name2taxid or scientific_name not in name2taxid:
            return {"kingdom": "Unknown", "phylum": "Unknown", "class": "Unknown"}
        taxid = name2taxid[scientific_name][0]
        lineage = ncbi.get_lineage(taxid)
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)
        taxonomy = {"kingdom": "Unknown", "phylum": "Unknown", "class": "Unknown"}
        for tid in lineage:
            rank = ranks.get(tid, "")
            if rank in taxonomy:
                taxonomy[rank] = names[tid]
        return taxonomy
    except Exception:
        return {"kingdom": "Unknown", "phylum": "Unknown", "class": "Unknown"}

# ---------------- Diversity Indices ----------------
def shannon_index(counts):
    """Shannon Diversity Index."""
    total = sum(counts.values())
    return -sum((c / total) * math.log(c / total) for c in counts.values() if c > 0)

def simpson_index(counts):
    """Simpson Diversity Index."""
    total = sum(counts.values())
    return 1 - sum((c / total) ** 2 for c in counts.values() if c > 0)

def chao1(counts):
    """Chao1 species richness estimator."""
    counts_list = list(counts.values())
    S_obs = sum(1 for c in counts_list if c > 0)
    F1 = sum(1 for c in counts_list if c == 1)
    F2 = sum(1 for c in counts_list if c == 2)
    if F2 == 0:
        return S_obs + (F1 * (F1 - 1)) / 2
    return S_obs + (F1 ** 2) / (2 * F2)
