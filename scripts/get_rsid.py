import sys

def get_dbsnp(dbsnp_file):
        
        """input: dbsnp_file chr, pos, alt, ref, rsid
        output: a dictonary of snp, key:chr_pos_alt_ref, value:rsid"""
        
        snps = {}
        f = open(dbsnp_file)
        for line in f:
                items = line.split()
                if len(items) ==5:
                        snp_id = "_".join([items[i] for i in [0, 1, 3, 4]])
                        snps[snp_id] = items[2]
                else:
                        continue
        return snps


def match_snp(snp_info, snps):
        """input: a list of snps[chr, pos, ref, alt], snp_dictionary
        output: rsid"""
        name1 = "_".join(snp_info)
        name2 = "_".join([snp_info[0], snp_info[1],snp_info[3], snp_info[2]])
        if name1 in snps:
                return snps[name1]
        
        elif name2 in snps:
                print(names)
                return snps[name2]
        
        else:
                return "NA"
       


def write_output(annotation, snps, output):
        """input annotation file, snp dictionary, output file name
        output a file the last row is filled with rsid"""

        out = open(out_anno, "w")
        out.write("\t".join(["chr","pos","varID","ref_vcf","alt_vcf","rsid"])+"\n") 
        f_anno = open(annotation)
        for line in f_anno:
                items = line.split()
               # if len(items) == 8:
               #         out.write(line)
               # else:
                        #print(line)
                info = [items[i] for i in [0, 1, 3, 4]]      
                rsid = match_snp(info, snps)
                #print([items[i] for i in [0, 1, 3, 4]])
                items[5]=rsid
                out.write("\t".join(items[:6])+"\n")
        out.close()




if __name__ == '__main__':
        dbsnp_file = sys.argv[1]
        annot_file = sys.argv[2]
        out_dir = sys.argv[3]
        out_prefix = sys.argv[4]
        out_anno = out_dir + out_prefix +".txt"
        print("getting snp dictionary")
        snps = get_dbsnp(dbsnp_file)
        print(len(snps))
        #print(snps.keys()[0])
        #print(match_snp(["1", "88338", "G", "A"], snps))
        print("writting new output")
        write_output(annot_file, snps, out_anno)
