from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction
import re

def fetch_gene_from_ncbi(ncbi_id, user_email):
    """Fetches a genetic sequence from NCBI using Entrez."""
    print(f"Fetching gene ID: {ncbi_id} from NCBI...")
    Entrez.email = user_email
    handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record.seq

def find_crispr_targets(dna_sequence):
    """Finds 20-bp CRISPR-Cas9 targets upstream of an 'NGG' PAM sequence."""
    print("Scanning DNA for CRISPR PAM (NGG) targets...")
    dna_string = str(dna_sequence)
    pam_pattern = re.compile(r'(?=([ATCG]GG))')
    
    guides =[]
    
    for match in pam_pattern.finditer(dna_string):
        start = match.start()
        # Ensure there is enough sequence upstream to grab a 20-bp guide
        if start >= 20:
            guide_seq = dna_string[start-20 : start]
            pam_seq = dna_string[start : start+3]
            
            # Use Biopython to calculate GC%
            gc_percent = gc_fraction(guide_seq) * 100
            
            # Biological constraint: We only want guides with 40-60% GC content
            if 40 <= gc_percent <= 60:
                guides.append({
                    'location': start - 20,
                    'guide_rna': guide_seq,
                    'pam': pam_seq,
                    'gc_content': round(gc_percent, 2)
                })
                
    return guides

if __name__ == "__main__":
    # 1. Fetch Human HBB (Hemoglobin Beta) - NM_000518
    # Replace with your actual email!
    my_email = "your_email@example.com" 
    gene_seq = fetch_gene_from_ncbi("NM_000518", my_email)
    
    # 2. Find and Filter Guide RNAs
    optimal_guides = find_crispr_targets(gene_seq)
    
    # 3. Output the Top 5 most efficient Guides
    print(f"\nSuccess! Found {len(optimal_guides)} OPTIMAL guide RNAs (40-60% GC Content).")
    print("--- Top 5 Guide RNAs ---")
    for i, guide in enumerate(optimal_guides[:5]):
        print(f"Option {i+1}:")
        print(f"   Guide: {guide['guide_rna']} | PAM: {guide['pam']}")
        print(f"   GC %:  {guide['gc_content']}% | Pos: {guide['location']}")
        print("-" * 30)