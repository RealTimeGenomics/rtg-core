/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.util;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;

import com.rtg.util.Counter;
import com.rtg.util.MultiSet;
import com.rtg.util.Pair;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;

/**
 * Accumulates proportions of child GT calls for each mother/father GT combination
 */
public final class GenotypeProportions {

  /**
   * Constructor
   */
  public GenotypeProportions() { }

  private final TreeMap<Parents, MultiSet<Genotype>> mAlleleSets = new TreeMap<>(Parents.PARENTS_COMPARATOR);


  /**
   * @param parent1 genotype for first parent
   * @param parent2 genotype for second parent
   * @param child genotype for child
   */
  public void addRecord(Genotype parent1, Genotype parent2, Genotype child) {
    final Parents parents = new Parents(parent1, parent2);
    if (!mAlleleSets.containsKey(parents)) {
      mAlleleSets.put(parents, new MultiSet<>(new TreeMap<Genotype, Counter>(Genotype.GENOTYPE_COMPARATOR)));
    }
    final MultiSet<Genotype> childAlleles = mAlleleSets.get(parents);
    childAlleles.add(child);
  }

  /**
   * Writes the output to the given stream
   * @param app destination for results
   * @throws IOException if an IO error occurs
   */
  public void writeResults(Appendable app) throws IOException {
    final List<Pair<Long, StringBuilder>> results = new ArrayList<>();
    for (Map.Entry<Parents, MultiSet<Genotype>> e : mAlleleSets.entrySet()) {
      final StringBuilder res = new StringBuilder();
      final Parents p = e.getKey();
      res.append("Parent GT:");
      for (Genotype i : p.mParents) {
        res.append(" ");
        res.append(i);
      }

      res.append(StringUtils.LS);
      res.append("Children:").append(StringUtils.LS);
      final MultiSet<Genotype> counts = e.getValue();
      final long tot = counts.totalCount();
      for (Genotype i : counts.keySet()) {
        res.append(String.format(Locale.ROOT, "%6s %10d %6.2f%%", i.toString(), counts.get(i), (double) counts.get(i) * 100.0d / tot));
        res.append(StringUtils.LS);
      }
      res.append("Total: ").append(Long.toString(tot)).append(StringUtils.LS);
      res.append(StringUtils.LS);
      results.add(new Pair<>(tot, res));
    }
    Collections.sort(results, new CountPairComparator());
    for (Pair<Long, StringBuilder> p : results) {
      app.append(p.getB().toString());
    }
  }

  //Inner Classes
  private static class Parents {
    static final ParentsComparator PARENTS_COMPARATOR = new Parents.ParentsComparator();

    Genotype[] mParents;
    Parents(Genotype parent1, Genotype parent2) {
      mParents = new Genotype[] {parent1, parent2};
      Arrays.sort(mParents, new Genotype.GenotypeComparator());
    }

    public boolean equals(Object o) {
      if (o == null || !(o instanceof  Parents)) {
        return false;
      }
      return Arrays.equals(mParents, ((Parents) o).mParents);
    }

    public int hashCode() {
      int hash = 0;
      for (Genotype i : mParents) {
        hash = Utils.pairHash(hash, i.hashCode());
      }
      return hash;
    }

    //If this modified to be non-thread safe, then make sure to remove the static instance
    private static class ParentsComparator implements Comparator<Parents>, Serializable {
      @Override
      public int compare(Parents a, Parents b) {
        int alleleCountA = 0;
        for (Genotype i : a.mParents) {
          alleleCountA += i.length();
        }
        int alleleCountB = 0;
        for (Genotype i : b.mParents) {
          alleleCountB += i.length();
        }
        if (alleleCountA < alleleCountB) {
          return 1;
        } else if (alleleCountA > alleleCountB) {
          return -1;
        }
        for (int i = 0; i < a.mParents.length && i < b.mParents.length; i++) {
          if (Genotype.GENOTYPE_COMPARATOR.compare(a.mParents[i], b.mParents[i]) < 0) {
            return -1;
          } else if (Genotype.GENOTYPE_COMPARATOR.compare(a.mParents[i], b.mParents[i]) > 0) {
            return 1;
          }
        }
        if (a.mParents.length < b.mParents.length) {
          return -1;
        } else if (a.mParents.length > b.mParents.length) {
          return 1;
        }
        return 0;
      }
    }
  }

  private static class CountPairComparator implements Comparator<Pair<Long, StringBuilder>>, Serializable {
    @Override
    public int compare(Pair<Long, StringBuilder> a, Pair<Long, StringBuilder> b) {
      return b.getA().compareTo(a.getA());
    }
  }


  private static void useReader(File f) throws IOException {
    final GenotypeProportions prop = new GenotypeProportions();
    try (VcfReader vcfR = VcfReader.openVcfReader(f)) {
      while (vcfR.hasNext()) {
        final VcfRecord rec = vcfR.next();
        final ArrayList<String> sampleGts = rec.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE);
        prop.addRecord(new Genotype(sampleGts.get(0)), new Genotype(sampleGts.get(1)), new Genotype(sampleGts.get(2)));
      }
    }
    prop.writeResults(System.err);
  }

  /**
   * @param args should have length 1 and be the filename of the file to process
   * @throws Exception could be anything
   */
  public static void main(String[] args) throws Exception {
    useReader(new File(args[0]));
  }

}
