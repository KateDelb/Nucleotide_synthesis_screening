#!/usr/bin/env python3
"""
HMMER domain screening for translated protein sequences
Identifies conserved functional domains that might indicate virulence/danger
"""

import pyhmmer
from pyhmmer.plan7 import HMM
import tempfile

def translate_sequence(dna_seq: str) -> str:
    """Translate DNA to protein using simple codon table"""
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    
    dna_seq = dna_seq.upper()
    protein = []
    
    # Start from first ATG
    start = dna_seq.find('ATG')
    if start == -1:
        start = 0
    
    for i in range(start, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        if len(codon) == 3:
            protein.append(codon_table.get(codon, 'X'))
    
    return ''.join(protein)

def run_hmmscan(protein_seq: str, db_path: str = "Pfam-A.hmm") -> dict:
    """
    Run HMMER hmmscan on protein sequence using pyhmmer
    
    Args:
        protein_seq: Protein sequence (string)
        db_path: Path to Pfam HMM database file
    
    Returns:
        dict with domain hits and risk assessment
    """
    
    try:
        from pathlib import Path
        
        db_file = Path(db_path)
        if not db_file.exists():
            return {
                "status": "DB_NOT_FOUND",
                "message": f"Pfam database not found at {db_path}. Download from: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
                "domains": [],
                "risk_flags": ["⚠️  Pfam database not found"]
            }
        
        # Load all HMMs from database into memory
        hmms = []
        with pyhmmer.plan7.HMMFile(str(db_file)) as hmm_file:
            for hmm in hmm_file:
                hmms.append(hmm)
        
        if not hmms:
            return {
                "status": "ERROR",
                "message": "No HMMs found in database file",
                "domains": [],
                "risk_flags": []
            }
        
        # Get alphabet from first HMM
        alphabet = hmms[0].alphabet
        
        # Create digital sequence properly
        text_seq = pyhmmer.easel.TextSequence(name=b"query", sequence=protein_seq.encode())
        digital_seq = text_seq.digitize(alphabet)
        
        domains = []
        
        # Search against all HMMs
        for hits in pyhmmer.hmmsearch(hmms, [digital_seq]):
            for hit in hits:
                if hit.included:  # Only keep included (significant) hits
                    domain = {
                        "name": hit.name.decode() if isinstance(hit.name, bytes) else hit.name,
                        "accession": hit.accession.decode() if isinstance(hit.accession, bytes) else hit.accession,
                        "evalue": hit.evalue,
                        "score": hit.score,
                        "description": hit.description.decode() if isinstance(hit.description, bytes) else hit.description
                    }
                    domains.append(domain)
        
        # Sort by E-value (best matches first)
        domains.sort(key=lambda x: x['evalue'])
        
        risk_flags = assess_domain_risk(domains)
        
        return {
            "status": "SUCCESS",
            "domains": domains,
            "risk_flags": risk_flags,
            "total_hits": len(domains)
        }
    
    except ImportError:
        return {
            "status": "PYHMMER_NOT_INSTALLED",
            "message": "pyhmmer not found. Install with: pip install pyhmmer",
            "domains": [],
            "risk_flags": ["⚠️  Cannot perform domain screening without pyhmmer"]
        }
    except Exception as e:
        return {
            "status": "ERROR",
            "message": str(e),
            "domains": [],
            "risk_flags": []
        }

def assess_domain_risk(domains: list) -> list:
    """
    Assess risk based on detected domains
    Look for dangerous Pfam domains
    """
    
    # Dangerous domain keywords
    danger_keywords = {
        "polymerase": "RNA/DNA polymerase — viral replication machinery",
        "reverse_transcriptase": "Reverse transcriptase — retroviral machinery",
        "protease": "Protease — viral/bacterial protein processing",
        "toxin": "Toxin domain — direct harm potential",
        "hemolytic": "Hemolytic toxin — destroys red blood cells",
        "botulinum": "Botulinum toxin — neurotoxin",
        "ricin": "Ricin — protein synthesis inhibitor",
        "shiga": "Shiga toxin — pathogenic",
        "anthrax": "Anthrax related — weaponizable",
        "integrase": "Integrase — retroviral integration",
        "envelope": "Envelope protein — viral surface protein",
        "spike": "Spike protein — viral attachment",
        "capsid": "Capsid — viral structural protein",
        "nucleocapsid": "Nucleocapsid — viral structural protein",
        "kinase": "Kinase — cell signaling (context-dependent)",
        "adhesin": "Adhesin — bacterial attachment/virulence",
        "invasin": "Invasin — bacterial invasion",
        "immune_evasion": "Immune evasion — evades host immunity",
    }
    
    risk_flags = []
    
    for domain in domains:
        domain_name = domain['name'].lower()
        score = domain['score']
        evalue = domain['evalue']
        
        # Check against danger keywords
        for keyword, description in danger_keywords.items():
            if keyword in domain_name:
                if evalue < 0.01 and score > 20:  # High confidence hit
                    risk_flags.append(f"🔴 HIGH: {description} (Pfam: {domain['name']}, E-value: {evalue:.2e})")
                elif evalue < 0.1 and score > 15:  # Medium confidence
                    risk_flags.append(f"🟡 MEDIUM: {description} (Pfam: {domain['name']}, E-value: {evalue:.2e})")
                break
    
    return risk_flags

def screen_with_hmmer(dna_seq: str) -> dict:
    """
    Complete screening: translate → run HMMER → assess risk
    """
    
    # Translate
    protein = translate_sequence(dna_seq)
    
    if len(protein) < 10:
        return {
            "error": "Protein too short (<10 aa) for domain screening",
            "protein_length": len(protein),
            "domains": [],
            "risk_flags": []
        }
    
    # Run HMMER
    hmmer_result = run_hmmscan(protein)
    
    return {
        "protein_length": len(protein),
        "protein_seq": protein,
        "hmmer_status": hmmer_result.get("status"),
        "domains_found": hmmer_result.get("total_hits", 0),
        "domains": hmmer_result.get("domains", []),
        "risk_flags": hmmer_result.get("risk_flags", []),
        "message": hmmer_result.get("message", "")
    }

if __name__ == "__main__":
    # Test with mystery sequence
    mystery = "ATGAATCAACATAATACCCAAATCAATAAATTTATCTTTTTTAGTAGCTTCAGAATCATCAGTACCACCCCAATTCAATTTATCAACAATAATAGGACCAGCAACAGCAATAGCAACATATTCATCATTAGAAGCAGGAGAATTCAAAACAATACCAACAGCTTTTTTATTATCAGCAGTATTCCAAGTTTTACAAACAGCACCATTAGAATTACCTTTTTCACAATCACCAGCATGCAAATTATCAATACCAGAATTTTTATGAGATCTAATATGAACCAAACCCAT"
    
    print("\n" + "="*70)
    print("HMMER DOMAIN SCREENING")
    print("="*70)
    
    result = screen_with_hmmer(mystery)
    
    print(f"\nProtein length: {result['protein_length']} aa")
    print(f"Protein: {result['protein_seq'][:50]}...")
    
    print(f"\nHMMER Status: {result['hmmer_status']}")
    
    if result['message']:
        print(f"Message: {result['message']}")
    
    print(f"\nDomains found: {result['domains_found']}")
    
    if result['domains']:
        print("\nDetailed Domain Hits:")
        for domain in result['domains'][:5]:
            print(f"  • {domain['name']}")
            print(f"    └─ Accession: {domain['accession']}")
            print(f"    └─ E-value: {domain['evalue']:.2e} | Score: {domain['score']:.1f}")
            if 'description' in domain:
                print(f"    └─ {domain['description']}")
    else:
        print("  No significant domain matches found")
    
    if result['risk_flags']:
        print("\nRisk Flags:")
        for flag in result['risk_flags']:
            print(f"  {flag}")
    else:
        print("\n  ✓ No dangerous domains detected")
    
    print("\n" + "="*70)
