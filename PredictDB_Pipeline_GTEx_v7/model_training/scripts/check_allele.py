import itertools 
ann1="/projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/input/annotations/snp_annotation/annotation4.txt"
ann2="/projects/b1047/zhong/Hepatocyte_project/data/genotype/./baseline.genotype.eqtl.excludeoutlier.header.vcf.anno"
with open(ann1) as f1, open(ann2) as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()
        print(len(lines1))
        print(len(lines2))
        i = 0
        j = 0
        while i < len(lines1) and j <len(lines2):
                while lines1[i][0]=='#' or lines1[i][0]=="c":
                        i += 1
                while lines2[j][0]== '#':
                        j += 1
                if lines1[0] != '#' and lines2[0] != '#':
                        items1 = lines1[i].split()
                        items2 = lines2[j].split()
                        
                        chr1,pos1,ref1,alt1 = int(items1[0]), int(items1[1]), items1[3], items1[4]
                        chr2, pos2,ref2,alt2 =  int(items2[0][3:]), int(items2[1]), items2[3], items2[4]
                        while chr1 < chr2:
                                i+=1
                        while chr1 > chr2:
                                j+=1
                        if chr1 == chr2:
                                while pos1 < pos2:
                                        i += 1
                                while pos1 > pos2:
                                        j += 1
                                if pos1 == pos2:
                                        if ref1==ref2 and alt1 == alt2:
                                                i+=1
                                                j+=1
                                        else:
                                                print(items1, items2)
                                                i+=1
                                                j+=1

                        
print(i,j)
                        






