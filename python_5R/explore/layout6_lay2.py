# 注：此脚本在本次会话中从 scratchpad 运行，依赖未入库的中间文件，路径需自行替换；结果见 docs/primer_design/layout6_*
import json, sys, numpy as np, pandas as pd
sys.path.insert(0,'/home/user/16S_5R/python_5R/explore')
import primer_design as pdz, compact_primers as cp
SC='/tmp/claude-0/-home-user-16S-5R/0a157df1-6424-5d79-bd19-c9f07612faab/scratchpad/'
MAIN={'A1-F':[8,['AGRGTTTGATYMTGGCTCA','YGWGTTTGATCCTGGCKCA','MGWGTTTGATCCTGGCTKA']],
'A1-R':[246,['TTACCYCACCAACWARCT','TTACCCCRCCAACTABCT','CGTTACCYTACCAACTARYT','CATTACCCYACCAACTARYT']],
'A2-F':[314,['CAYKGGGACTGAGACACKG','CACTGGRACTGARAYACGG','CARKGGRACTGAGACACGG']],
'A2-R':[504,['GCTGGCACGKARTTAGCCR','GCTGGCACGAASTTAGYCS','GCTGGCACAKAGTTAGYYG']],
'A3-F':[518,['CCAGCAGCCGCGGTAABWC','CCAGCAGCCGCGGTAAAAC','CCASCASCCGCRGTAATAC']],
'A3-R':[774,['TCTAATCCYGTTYGCTMCC','TCTRATCGTCTTCGAWCCY','TCTAATCCKGTTYGMTCCC']],
'A4-F':[788,['TTAGAKACCCYGGTAGTCY','TYAGATACCGTCGTAGTYY','TTAGAWACCCSRGTAGTCC']],
'A4-R':[954,['GCGTWGCATCGAATTAADC','GCGTWKCKTCGAATTAAAC','GTGTWGCVTCGAATTAAAC']],
'A5-F':[967,['MWACGCGARGAACCTTACC','MWACMCGAAGAACCTTACC','RTACMCGAARAACCTTACC']],
'A5-R':[1174,['TCATCCCCACCTTCCTCCN','TCGTCCCCACCTTCCTCYD','TCATCCYCACCTTCCTCYR']],
'A6-F':[1222,['GGCTRCACACRTGCTACAA','GGCKACACACGTGMTACAA','GGCTWCACGCGTGCTACAW']],
'A6-R':[1486,['CTTGTTACGACTTCRYCCY','CTTGTTACGACTTYTMCWT','CTTGTTACGACTTMACHCC']]}
SUP={'A1-R':['GTTACCCCKCCAACWAGCT','GTTRCCCYGCCAACTAGCT','GTTACCCCGCCRWCWAGCT','GTTACCCCGCCWRCTRCCT','GTCACCCCGCCAWCAAGCT','ATTACCCCRCCGTCAAGCT'],
'A4-F':['TTAGATACCCTGGTAGTHC','TTAGATACCYTGGTAGTSC'],'A4-R':['GTGTATYRTCGAATTAARC','GCGTABCTTCRAATTAAAC'],
'A5-F':['CRACGCGAAGAACCTTAYC','CAACGCGAAGAACCTTWSC'],'A6-F':['GGCYWCACACGTACTACAA','GGCTTCACACGTCATACAA']}
oligos=[]
for s,(st,ps) in MAIN.items():
    for i,p in enumerate(ps,1): oligos.append(dict(name=f'16S-{s}.{i}',site=s,start=st,seq=p,kind='main'))
    for j,p in enumerate(SUP.get(s,[]),len(ps)+1): oligos.append(dict(name=f'16S-{s}.{j}s',site=s,start=st,seq=p,kind='supp'))
df=pd.DataFrame(oligos); df['nt']=df.seq.str.len(); df['exp']=df.seq.apply(lambda x:len(pdz.expand(x)))
tm=df.seq.apply(lambda x:[cp.oligo_tm(e)[0] for e in pdz.expand(x)]); df['tm_min']=tm.apply(min).round(1); df['tm_max']=tm.apply(max).round(1)
df['gc']=df.seq.apply(pdz.gc_content).round(0)
print(df[['name','seq','nt','exp','tm_min','tm_max']].to_string(index=False)); print(len(df),'oligos',df.exp.sum(),'expansions')
df.to_csv(SC+'lay2_oligos.csv',index=False)
# 池内二聚体（全部 oligo 同管）
by={r['name']:[r['seq']] for r in oligos}
sev,hp=cp.pool_issues(by)
print('严重二聚体',len(sev),'发夹',len(hp))
for x in sorted(sev,key=lambda t:t[4])[:30]: print(x)
for x in hp[:10]: print('hp',x)
