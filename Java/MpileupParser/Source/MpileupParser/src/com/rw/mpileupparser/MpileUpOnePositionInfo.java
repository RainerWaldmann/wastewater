/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package com.rw.mpileupparser;

import com.rw.mpileupparser.modifications.ModificationsBase;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

/**
 *
 * @author raine
 */
public class MpileUpOnePositionInfo {

    public final char refBase;
    private int matchesRefForward = 0;
    private int matchesRefReverse = 0;
    private int fwdQVsum = 0;
    private int revQVsum = 0;
    private int totalReadsOnPos = 0;
    Map<ModificationsBase, ModificationsBase> modifications = null;

    public int getRefDepth() {
        return matchesRefForward + matchesRefReverse;
    }

    public int getMatchesRefForward() {
        return matchesRefForward;
    }

    public int getMatchesRefReverse() {
        return matchesRefReverse;
    }

    public Optional<Float> getFWDqv() {
        return matchesRefForward == 0 ? Optional.empty() : Optional.of((float) fwdQVsum / matchesRefForward);
    }

    public Optional<Float> getREVqv() {
        return matchesRefReverse == 0 ? Optional.empty() : Optional.of((float) revQVsum / matchesRefReverse);
    }

    public int getTotalReadsOnPos() {
        return totalReadsOnPos;
    }
    
    public float getRefQV() {
        float sumQV = 0;
        int n = 0;
        Optional<Float> qV = getFWDqv();
        if (qV.isPresent()) {
            n++;
            sumQV += qV.get();
        }
        qV = getREVqv();
        if (qV.isPresent()) {
            n++;
            sumQV += qV.get();
        }
        return sumQV / n;
    }

    public MpileUpOnePositionInfo(char refBase) {
        this.refBase = refBase;
    }
   
    /**
     * if # or  totalreadsonpos must be incremented but refbases not -> use totalReadsOnpos to calc frequencies
     */
    public void incrementTotalReadsOnPos(){
        totalReadsOnPos++;
    }

    public void addRefMatch(boolean reverse, char qv) {
        if (reverse) {
            matchesRefReverse++;
            revQVsum += ModificationsBase.getQV.apply(qv);
        } else {
            matchesRefForward++;
            fwdQVsum += ModificationsBase.getQV.apply(qv);
        }
    }

    public void addModification(ModificationsBase mod, boolean reverse, char qv) {
        if (modifications == null) {
            modifications = new HashMap<>();
            modifications.put(mod, mod);
        } else {
            ModificationsBase m = modifications.get(mod);
            if (m == null) {
                modifications.put(mod, mod);
            } else {
                m.incrementCounts(reverse, qv);
            }
        }
    }

}
