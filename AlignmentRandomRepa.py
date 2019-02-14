from AlignmentRepa import *

# historyRepasShuffle_u :: HistoryRepa -> Int -> HistoryRepa

def historyRepasShuffle_u(aa,s):
    (vaa,maa,saa,raa) = aa
    (n,_) = raa.shape
    np.random.seed(s)
    ll = []
    for i in range(n):
        aa = np.copy(raa[i])
        np.random.shuffle(aa)
        ll.append(aa)
    rbb = np.stack(ll)
    return (vaa,maa,saa,rbb)

# systemsDecompFudsHistoryRepasMultiplyWithShuffle :: 
#   Integer -> Integer -> System -> DecompFud -> HistoryRepa -> Map (State,Fud) (HistoryRepa,HistoryRepa)

def systemsDecompFudsHistoryRepasMultiplyWithShuffle(mult,seed,uu,df,aa):
    def unit(ss):
        return setStatesHistogramUnit(sset([ss]))
    aahh = histogramsHistory
    fder = fudsDerived
    empty = historyRepaEmpty
    size = historyRepasSize
    vars = historyRepasSetVariable
    hhhr = systemsHistoriesHistoryRepa
    def select(uu,ss,hh):
        return historyRepasHistoryRepasHistoryRepaSelection_u(hhhr(uu,aahh(unit(ss))),hh)
    def hrhrred(hr,vv):
        return setVarsHistoryRepasHistoryRepaReduced(vv,hr)
    applyFud = systemsFudsHistoryRepasMultiply_u
    hrconcat = vectorHistoryRepasConcat_u
    hrshuffle = historyRepasShuffle_u
    def shuffle(xx):
        z = size(xx)
        if z == 0:
            return empty()
        return hrconcat([hrshuffle(xx,seed+i*z) for i in range(1,mult+1)])
    def apply(zz,vv,aa):
        qq = sdict()
        for ((ss,ff),yy) in zz.items():
            aa1 = select(uu,ss,aa)
            aaxx1 = shuffle(aa1)
            ww = fder(ff)
            bb = empty()
            if size(aa1) > 0:
                bb = hrhrred(applyFud(uu,ff,aa1),vv|ww)
            if size(aaxx1) > 0:
                bbxx = hrhrred(applyFud(uu,ff,aaxx1),vv|ww)
            qq[(ss,ff)] = (bb,bbxx)
            qq.update(apply(yy,vv,bb))
        return qq
    return apply(df,vars(aa),aa)

