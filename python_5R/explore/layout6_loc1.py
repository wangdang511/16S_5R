# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/layout6_*
import pickle, numpy as np, pandas as pd, warnings, collections, sys
warnings.filterwarnings('ignore'); pd.set_option('display.width',250); pd.set_option('display.max_rows',300)
import primer_design as pdz
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
B=pickle.load(open(SC+'big_labeled.pkl','rb')); S=pickle.load(open(SC+'S.pkl','rb'))
ids=B['ids']; tax=B['tax']
ph=np.array([tax[i][1] if tax[i][1] else 'unclassified' for i in ids])
G=pdz.Reference(B['aln'],ids,np.array(['Bacteria']*len(ids)),ph,S.col_of_pos,S.ref_seq,{})
genus=np.array([tax[i][5] if len(tax[i])>5 and tax[i][5] else '' for i in ids])
SITES={('F',8):['AGRGTTTGATYMTGGCTCA','YGWGTTTGATCCTGGCKCA','MGWGTTTGATCCTGGCTKA'],
('R',246):['TTACCYCACCAACWARCT','TTACCCCRCCAACTABCT','CGTTACCYTACCAACTARYT','CATTACCCYACCAACTARYT'],
('F',314):['CAYKGGGACTGAGACACKG','CACTGGRACTGARAYACGG','CARKGGRACTGAGACACGG'],
('R',504):['GCTGGCACGKARTTAGCCR','GCTGGCACGAASTTAGYCS','GCTGGCACAKAGTTAGYYG'],
('F',518):['CCAGCAGCCGCGGTAABWC','CCAGCAGCCGCGGTAAAAC','CCASCASCCGCRGTAATAC'],
('R',774):['TCTAATCCYGTTYGCTMCC','TCTRATCGTCTTCGAWCCY','TCTAATCCKGTTYGMTCCC'],
('F',788):['TTAGAKACCCYGGTAGTCY','TYAGATACCGTCGTAGTYY','TTAGAWACCCSRGTAGTCC'],
('R',954):['GCGTWGCATCGAATTAADC','GCGTWKCKTCGAATTAAAC','GTGTWGCVTCGAATTAAAC'],
('F',967):['MWACGCGARGAACCTTACC','MWACMCGAAGAACCTTACC','RTACMCGAARAACCTTACC'],
('R',1174):['TCATCCCCACCTTCCTCCN','TCGTCCCCACCTTCCTCYD','TCATCCYCACCTTCCTCYR'],
('F',1222):['GGCTRCACACRTGCTACAA','GGCKACACACGTGMTACAA','GGCTWCACGCGTGCTACAW'],
('R',1486):['CTTGTTACGACTTCRYCCY','CTTGTTACGACTTYTMCWT','CTTGTTACGACTTMACHCC']}
L=19
hit={}; dat={}; MAT={}
for (o,st),ps in SITES.items():
    h=np.zeros(len(ids),bool); Lm=max(len(p) for p in ps)
    M=pdz.site_matrix(G,st,L); MAT[(o,st)]=M; dat[(o,st)]=(M!=0).all(1)&(M!=ord('-')).all(1)
    for p in ps:
        Mp=pdz.site_matrix(G,st,len(p)); dd=(Mp!=0).all(1)&(Mp!=ord('-')).all(1)
        top=p if o=='F' else pdz.revcomp(p); side='right' if o=='F' else 'left'
        hh=np.zeros(len(ids),bool); hh[dd]=pdz._match(Mp[dd],top,side); h|=hh
    hit[(o,st)]=h
pickle.dump((hit,dat),open(SC+'loc_hit.pkl','wb'))
TARGETS={'Bifidobacterium':['F8','R246'],'Prevotella':['F8','R246'],'Roseburia':['F8','R246'],'Dialister':['F8','R246'],'Blautia':['F8','R246'],'Megasphaera':['F8','R246'],'Veillonella':['F8','R246','F314','R504','R954','F967','R1174','F788'],
'Mycoplasma':['R954','F788','F967','R1174','R246','F314'],'Ureaplasma':['R954','F788'],'Streptococcus':['R954','F967','R1174','F788','R504','F314'],'Staphylococcus':['R954','F788','R246','F967','R1174'],
'Neisseria':['F1222','R1486','F518'],'Oscillospira':['F1222','R1486'],'Faecalibacterium':['F1222','R1486','F314','R504'],'Ruminococcus':['F1222','R1486'],
'Haemophilus':['F314','R504','R246','F8'],'Peptoniphilus':['F314','R504'],'Enterococcus':['F314','R504','F518','F967','R1174'],'Gemella':['F314','R504','F788','R954','F967'],'Lactobacillus':['R954','F518','F788']}
print('各属在各位点的匹配率（%）；n=有数据序列数')
for g,sites in TARGETS.items():
    m=genus==g; line='%-16s n=%3d '%(g,m.sum())
    for s in sites:
        k=(s[0],int(s[1:])); mm=m&dat[k]
        line+=' %s:%s'%(s,('%d(%d)'%(round(100*hit[k][mm].mean()),mm.sum()) if mm.sum()>=3 else '–'))
    print(line)
print('\n失配序列的典型位点序列（19-mer，正链方向）与最近引物的错配位置')
def mm_pos(top,seq):
    return [i for i,(c,b) in enumerate(zip(top,seq)) if b not in pdz.IUPAC[c]]
for g,sites in TARGETS.items():
    m=genus==g
    for s in sites:
        k=(s[0],int(s[1:])); sel=m&dat[k]&~hit[k]
        if sel.sum()<3: continue
        seqs=[bytes(x).decode() for x in MAT[k][sel]]
        top=[p if k[0]=='F' else pdz.revcomp(p) for p in SITES[k]]
        c=collections.Counter(seqs).most_common(2)
        for q,nq in c:
            best=min(top,key=lambda t:len(mm_pos(t,q)) if len(t)==L else len(mm_pos(t,q[:len(t)])) if k[0]=='F' else len(mm_pos(t,q[-len(t):])))
            q2=q if len(best)==L else (q[:len(best)] if k[0]=='F' else q[-len(best):])
            print('%-16s %-6s 失配%d/%d  %s  错配位(0起,正链)%s'%(g,s,sel.sum(),(m&dat[k]).sum(),q,mm_pos(best,q2)))
