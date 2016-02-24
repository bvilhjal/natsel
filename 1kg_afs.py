import h5py
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt 
import scipy as sp

kg_dir = '/faststorage/project/TheHonestGene/1Kgenomes/phase3/'

def gen_anc_afs_plot(ancestry = 'AFR', plot_prefix = '/Users/bjv/Dropbox/Cloud_folder/tmp/1kg_AFS', data_filter=0.01):
    h5f = h5py.File('%s1k_genomes_hg.hdf5'%kg_dir)
    
    ancestries = sp.unique(h5f['indivs']['ancestry'][...])
    #ancestries = sp.unique(h5f['indivs']['continent'][...])
    
    for ancestry in ancestries:
        anc_filter = h5f['indivs']['ancestry'][...]==ancestry
        #anc_filter = h5f['indivs']['continent'][...]==ancestry
        
        acs = []
        for chrom in range(1,23):
            print 'Working on chromosome %d'%chrom
            chr_str = 'chr%d'%chrom
            snps  = sp.array(h5f[chr_str]['calldata']['snps'][...],dtype='int8')
            snps = snps[:,anc_filter]
            rand_filt = sp.random.random(len(snps))
            rand_filt = rand_filt<data_filter
            snps = snps[rand_filt]
            print 'SNPs loaded and filtered, leaving %d SNPs'%len(snps)
            (m,n) = snps.shape
            ac = sp.sum(snps,1)
            flip_filter = ac>n
            ac[flip_filter]=2*n-ac[flip_filter] 
            ac = ac[ac>0]
            ac = ac[ac<30]
            acs.extend(ac)
            
        acs = sp.array(acs,dtype='int')
        min_ac = acs.min()
        max_ac = acs.max()
        sp.bincount(acs)
        plt.clf()
        plt.hist(acs, bins=sp.arange(min_ac-0.5, max_ac + 1.5, 1))
        plt.title('%s AFS'%ancestry)
        plt.savefig('%s_%s.png'%(plot_prefix,ancestry))
    
    
if __name__=='__main__':
#     for ancestry in ['EUR','SAS','EAS','AFR','AMR']:
#         gen_anc_afs_plot(ancestry)
    for ancestry in ['EUR','SAS','EAS','AFR','AMR']:
        in_file = '/Users/bjarnivilhjalmsson/Dropbox/1KG_share_w_Asger/1kg_AFS_all_%s.txt'%ancestry
        out_file = '/Users/bjarnivilhjalmsson/Dropbox/1KG_share_w_Asger/1kg_AFS_all_%s_counts.txt'%ancestry
        with open(in_file) as f:
            l = (f.next()).split()
            n = int(l[1])
            count_dict = dict(zip(range(n+1),[0]*(n+1)))
            print f.next()
            for line in f:
                l = line.split()
                count_dict[int(l[1])]+=1
            with open(out_file,'w') as outf:
                outf.write('AC num_SNPs\n') 
                for i in range(n+1):
                    outf.write('%d %d\n'%(i,count_dict[i])) 
                
