import h5py
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt 
import scipy as sp
from itertools import izip

# kg_dir = '/faststorage/project/TheHonestGene/1Kgenomes/phase3/'
kg_dir = '/Users/bjv/Data/1Kgenomes/phase3/'

def gen_anc_afs_plot(ancestry = 'AFR', plot_prefix = '/Users/bjv/Dropbox/Cloud_folder/tmp/1kg_AFS_all', 
                     outfile_prefix='/Users/bjv/Dropbox/Cloud_folder/tmp/1kg_AFS_all',
                     data_filter=1, chunk_size = 10000):
    h5f = h5py.File('%s1k_genomes_hg.hdf5'%kg_dir)
    
    anc_filter = h5f['indivs']['continent'][...]==ancestry
    num_indivs = sp.sum(anc_filter)
    print "%d individuals with %s ancestry are used"%(num_indivs,ancestry)
    
    acs = []
    sids_list = [] 
    for chrom in range(1,23):
        print 'Working on chromosome %d'%chrom
        chr_str = 'chr%d'%chrom
        num_snps = len(h5f[chr_str]['calldata']['snps'])
        assert num_snps==len(h5f[chr_str]['variants']['ID']), 'WTF?'
        for start_i in range(0,num_snps,chunk_size):
            snps  = sp.array(h5f[chr_str]['calldata']['snps'][start_i:start_i+chunk_size],dtype='int8')
            sids  = h5f[chr_str]['variants']['ID'][start_i:start_i+chunk_size]
            snps = snps[:,anc_filter]
            if data_filter<1:
                rand_filt = sp.random.random(len(snps))
                rand_filt = rand_filt<data_filter
                snps = snps[rand_filt]
                sids = sids[rand_filt]
            (m,n) = snps.shape
            ac = sp.sum(snps,1)
            flip_filter = ac>n
            ac[flip_filter]=2*n-ac[flip_filter]
            #Plotting filter 
            acs.extend(ac)
            sids_list.extend(sids)
            if start_i%1000000==0:
                print 'Parsed %d SNPs'%start_i
        print '%d SNPs loaded and filtered'%num_snps
        
    h5f.close()
    
    print '%d ACs and SIDs found'%len(acs)

    print 'Storing the AFS'
    with open('%s_%s.txt'%(outfile_prefix,ancestry),'w') as f:
        f.write('# %d individuals used\n'%num_indivs)
        f.write('SID    AC\n')
        for sid, ac in izip(sids_list,acs):
            f.write('%s    %d\n'%(sid,ac))

    print 'Plot things'
    acs = sp.array(acs,dtype='int')
    acs = acs[acs>0]
    acs = acs[acs<30]
    min_ac = acs.min()
    max_ac = acs.max()
    sp.bincount(acs)
    plt.clf()
    plt.hist(acs, bins=sp.arange(min_ac-0.5, max_ac + 1.5, 1))
    plt.title('%s AFS'%ancestry)
    plt.savefig('%s_%s.png'%(plot_prefix,ancestry))
    
    
for ancestry in ['EUR','SAS','EAS','AFR','AMR']:
    gen_anc_afs_plot(ancestry)