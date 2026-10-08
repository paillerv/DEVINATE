import sys
import os
import pysam

def load_classification_status(classification_path):
    """
    Charge le dictionnaire des décisions de classification : {read_id: decision} (KEEP / DROP)
    """
    status_dict = {}
    if not classification_path or not os.path.exists(classification_path):
        return status_dict
        
    with open(classification_path, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 12:
                read_id = parts[0]
                decision = parts[11].strip()
                status_dict[read_id] = decision
    print(f"DEBUG: {len(status_dict)} reads chargés depuis le fichier de classification.")
    return status_dict

def load_insertions_mapping(insertions_path, classification_path=None):
    """
    Charge toutes les insertions individuelles, récupère leur statut (KEEP/DROP) 
    via classification.tsv, et les organise par chromosome.
    """
    chrom_to_insertions = {}
    class_dict = load_classification_status(classification_path)
    
    with open(insertions_path, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            if parts[0].lower() in ("chr", "chrom") or parts[1].lower() == "start":
                continue
                
            chrom = parts[0]
            if "," in chrom:
                chrom = chrom.split(",")[0]
                
            try:
                start = int(parts[1].split(",")[0]) if "," in parts[1] else int(parts[1])
                end = int(parts[2].split(",")[0]) if "," in parts[2] else int(parts[2])
            except ValueError:
                continue
                
            read_id = parts[3]
            decision = class_dict.get(read_id, "UNKNOWN")
            
            if chrom not in chrom_to_insertions:
                chrom_to_insertions[chrom] = []
            chrom_to_insertions[chrom].append((start, end, read_id, decision))
            
    return chrom_to_insertions, class_dict

def audit_locus_fixed_v12(bam_path, clusters_path, insertions_path, output_path, classification_path=None, window=500):
    chrom_to_insertions, class_dict = load_insertions_mapping(insertions_path, classification_path)
    bam = pysam.AlignmentFile(bam_path, "rb")
    results = []
    all_read_details = []  # Pour stocker les détails de chaque read pour l'export
    
    print("Audit VAF robuste avec export détaillé des reads en cours...")
    with open(clusters_path, 'r') as f:
        header = f.readline() # Ignore l'en-tête
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            confidence = parts[5]
            
            # 1. Récupération des reads KEEP validés par le fichier BED
            alt_reads_set = set() 
            keep_count = 0
            
            if chrom in chrom_to_insertions:
                for bed_start, bed_end, read_id, decision in chrom_to_insertions[chrom]:
                    if not (bed_end < start or bed_start > end):
                        if decision == "KEEP":
                            if read_id not in alt_reads_set:
                                alt_reads_set.add(read_id)
                                keep_count = keep_count + 1
                                all_read_details.append((chrom, start, end, read_id, "KEEP"))
            
            # 2. Analyse de la couverture locale dans le BAM (Attrape les DROP et les Ref)
            loc_start = max(0, start - window)
            loc_end = end + window
            
            drop_count = 0
            ref_clean = 0
            ref_suspicious = 0
            
            try:
                for read in bam.fetch(chrom, loc_start, loc_end):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    
                    r_id = read.query_name
                    
                    # Si c'est déjà un KEEP enregistré, on passe
                    if r_id in alt_reads_set and any(d[3] == r_id and d[4] == "KEEP" for d in all_read_details if d[0]==chrom and d[1]==start):
                        continue
                    
                    # Si le read du BAM est présent dans la classification avec le statut DROP
                    decision = class_dict.get(r_id, "UNKNOWN")
                    if decision == "DROP":
                        # Évite les doublons si le read est déjà compté pour ce locus
                        if not any(d[3] == r_id and d[4] == "DROP" for d in all_read_details if d[0]==chrom and d[1]==start):
                            drop_count += 1
                            alt_reads_set.add(r_id)
                            all_read_details.append((chrom, start, end, r_id, "DROP"))
                        continue
                    
                    # S'il est déjà dans alt_reads_set via une autre condition
                    if r_id in alt_reads_set:
                        continue
                    
                    # --- Vrais reads de référence potentiels ---
                    r_start = read.reference_start
                    r_end = read.reference_end
                    
                    spans_locus = (r_start is not None and r_end is not None and r_start <= start and r_end >= end)
                    
                    has_central_clip = False
                    cigar = read.cigartuples
                    if cigar:
                        curr_pos = r_start
                        for op, length in cigar:
                            if op in (0, 7, 8): # M, =, X
                                if curr_pos is not None:
                                    curr_pos += length
                            elif op == 4: # Soft clip
                                if curr_pos and abs(curr_pos - start) < 50:
                                    has_central_clip = True
                                    break
                    
                    if spans_locus and not has_central_clip:
                        ref_clean += 1
                        all_read_details.append((chrom, start, end, r_id, "ref_clean"))
                    else:
                        ref_suspicious += 1
                        all_read_details.append((chrom, start, end, r_id, "ref_suspicious"))
                        
            except ValueError:
                continue
           
            # Calculs globaux
            reads_classified = keep_count + drop_count
            total_ref = ref_clean + ref_suspicious
            total_coverage = reads_classified + total_ref
            
            # 3. Calculs des VAF
            ratio_suspicious = (ref_suspicious / total_ref) if total_ref > 0 else 0.0
            
            if reads_classified == 0 and ratio_suspicious > 0.70 and total_ref > 10:
                vaf_strict = 100.0
                vaf_with_drops = 100.0
            else:
                total_inf_strict = keep_count + ref_clean
                vaf_strict = (keep_count / total_inf_strict * 100) if total_inf_strict > 0 else 0.0
                
                total_inf_drops = keep_count + drop_count + ref_clean
                vaf_with_drops = ((keep_count + drop_count) / total_inf_drops * 100) if total_inf_drops > 0 else 0.0
                
            results.append({
                'chr': chrom, 'start': start, 'end': end,
                'confidence': confidence,
                'total_coverage': total_coverage,
                'reads_classified': reads_classified,
                'keep_count': keep_count,
                'drop_count': drop_count,
                'ref_clean': ref_clean,
                'ref_suspicious': ref_suspicious,
                'vaf_strict': round(vaf_strict, 2),
                'vaf_with_drops': round(vaf_with_drops, 2)
            })
            
    bam.close()
    
    # 1. Écriture du fichier principal de résumé par locus
    with open(output_path, 'w') as out:
        out.write("chr\tstart\tend\tconfidence\ttotal_coverage\treads_classified\tkeep_count\tdrop_count\tref_clean\tref_suspicious\tVAF_strict\tVAF_with_drops\n")
        for r in results:
            out.write(f"{r['chr']}\t{r['start']}\t{r['end']}\t{r['confidence']}\t"
                      f"{r['total_coverage']}\t{r['reads_classified']}\t{r['keep_count']}\t{r['drop_count']}\t"
                      f"{r['ref_clean']}\t{r['ref_suspicious']}\t"
                      f"{r['vaf_strict']}\t{r['vaf_with_drops']}\n")
                      
    # 2. Écriture du fichier secondaire avec le détail de chaque read
    reads_output_path = output_path.replace(".tsv", "_reads_detail.tsv")
    with open(reads_output_path, 'w') as out_reads:
        out_reads.write("chr\tstart\tend\tread_id\tcategory\n")
        for chrom, start, end, read_id, category in all_read_details:
            out_reads.write(f"{chrom}\t{start}\t{end}\t{read_id}\t{category}\n")
                      
    print(f"Terminé !\n - Rapport VAF : {output_path}\n - Détail des reads : {reads_output_path}")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python calculate_allele_frequency.py <global.bam> <clusters.tsv> <insertions.bed> <output.tsv> [classification.tsv]")
        sys.exit(1)
        
    bam_file = sys.argv[1]
    clusters_file = sys.argv[2]
    bed_file = sys.argv[3]
    output_file = sys.argv[4]
    class_file = sys.argv[5] if len(sys.argv) > 5 else None
    
    audit_locus_fixed_v12(bam_file, clusters_file, bed_file, output_file, class_file)
