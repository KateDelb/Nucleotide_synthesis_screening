
from hmmer_screening import screen_with_hmmer
import tempfile
import subprocess

MYSTERY_SEQ = "ATGAATCAACATAATACCCAAATCAATAAATTTATCTTTTTTAGTAGCTTCAGAATCATCAGTACCACCCCAATTCAATTTATCAACAATAATAGGACCAGCAACAGCAATAGCAACATATTCATCATTAGAAGCAGGAGAATTCAAAACAATACCAACAGCTTTTTTATTATCAGCAGTATTCCAAGTTTTACAAACAGCACCATTAGAATTACCTTTTTCACAATCACCAGCATGCAAATTATCAATACCAGAATTTTTATGAGATCTAATATGAACCAAACCCAT"

def validate_DNA(seq: str):
    valid = all(b in "AGCT" for b in seq.upper())
    if not valid:
        raise ValueError("Expected DNA sequence with only A G C T")
    return "Initial scan: sequence is valid."

def check_orf(seq: str):
    start_pos = None
    end_pos = None

    # Find first ATG
    for b in range(len(seq) - 2):
        if seq[b:b+3] == "ATG":
            start_pos = b
            break
    
    # If ATG found, look for stop codons in the same reading frame (triplets after ATG)
    if start_pos is not None:
        for b in range(start_pos + 3, len(seq) - 2, 3):  # Step by 3 (codons)
            if seq[b:b+3] in ["TAG", "TAA", "TGA"]:
                end_pos = b
                break

    return [start_pos, end_pos]

def classify_organism(seq: str):
    """Predict organism type from CpG ratio and GC content
    NOTE: this is an initial rough guess at the organism type and should not be seen as conclusive
    """
    seq = seq.upper()
    
    gc_content = (seq.count('G') + seq.count('C')) / len(seq)
    cpg_count = seq.count('CG')
    expected_cpg = (seq.count('C') / len(seq)) * (seq.count('G') / len(seq)) * len(seq)
    cpg_ratio = cpg_count / expected_cpg if expected_cpg > 0 else 0
    
    if cpg_ratio < 0.5:
        organism = "EUKARYOTIC/VIRAL (CpG depleted)"
    elif cpg_ratio > 0.8:
        organism = "BACTERIAL (normal CpG)"
    else:
        organism = "AMBIGUOUS"
    
    return {
        "gc_content": f"{gc_content * 100:.1f}%",
        "cpg_ratio": f"{cpg_ratio:.2f}",
        "prediction": organism
    }

def run_blastn(seq_str, db="nt"):
    """BLAST nucleotide search (checks for exact/near matches to known sequences)"""
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">query\n")
            f.write(seq_str + "\n")
            f.flush()
            
            result = subprocess.run(
                ["blastn", "-query", f.name, "-db", db, "-outfmt", "6", "-max_target_seqs", "5"],
                capture_output=True,
                text=True,
                timeout=30,
                check=False
            )
            
            if result.returncode == 0 and result.stdout.strip():
                hits = []
                for line in result.stdout.strip().split('\n'):
                    parts = line.split('\t')
                    hits.append({
                        "subject": parts[1],
                        "identity": f"{float(parts[2]):.1f}%",
                        "length": parts[3],
                        "evalue": parts[10]
                    })
                return {"status": "HITS_FOUND", "results": hits}
            else:
                return {"status": "NO_HITS", "message": "No significant matches in NCBI nt database"}
    except FileNotFoundError:
        return {"status": "BLASTN_NOT_INSTALLED", "message": "blastn not found. Install NCBI BLAST tools."}
    except Exception as e:
        return {"status": "ERROR", "message": str(e)}
 
def run_blastx(seq_str, db="nr"):
    """BLASTX protein search (translate sequence and search protein database)"""
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">query\n")
            f.write(seq_str + "\n")
            f.flush()
            
            result = subprocess.run(
                ["blastx", "-query", f.name, "-db", db, "-outfmt", "6", "-max_target_seqs", "5"],
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode == 0 and result.stdout.strip():
                hits = []
                for line in result.stdout.strip().split('\n'):
                    parts = line.split('\t')
                    hits.append({
                        "subject": parts[1],
                        "identity": f"{float(parts[2]):.1f}%",
                        "length": parts[3],
                        "evalue": parts[10]
                    })
                return {"status": "HITS_FOUND", "results": hits}
            else:
                return {"status": "NO_HITS", "message": "No significant matches in NCBI nr database"}
    except FileNotFoundError:
        return {"status": "BLASTX_NOT_INSTALLED", "message": "blastx not found. Install NCBI BLAST tools."}
    except Exception as e:
        return {"status": "ERROR", "message": str(e)}
 

if __name__ == "__main__":
    flags = []

    #validate seq
    validate_DNA(MYSTERY_SEQ)

    #Check for start and stop codons
    if check_orf(MYSTERY_SEQ)[0] is not None and check_orf(MYSTERY_SEQ)[1] is None: 
        flags.append("Truncated ORF, possibly trying to avoid detection")

    #Check for organism based on codon bias
    organism = classify_organism(MYSTERY_SEQ)
    if "VIRAL" in organism or "BACTERIAL" in organism:
        flags.append(f"Organism: {organism}")

    #Check for nucleotide similarities
    if run_blastn(MYSTERY_SEQ)['status'] == "NO_HITS":
        flags.append("No significant matches in NCBI nt database")

    #Check for amino acid similarities
    if run_blastn(MYSTERY_SEQ)['status'] == "NO_HITS":
            flags.append("No significant matches in NCBI nr (non-redundant protein) database")

    #Check for conserved dangerous domains in protein structure with hmmer
    hmmer_result = screen_with_hmmer(MYSTERY_SEQ)
    print(f"\n2. Domain Screening (HMMER):")
    print(f"   Status: {hmmer_result['hmmer_status']}")
    print(f"   Protein length: {hmmer_result['protein_length']} aa")
    print(f"   Domains found: {hmmer_result['domains_found']}")
    if hmmer_result["domains_found"] > 0:
        flags.append("Potential dangerous protein domains found with HMMer")

    print(f"{len(flags)} flags: {flags}")

    
        