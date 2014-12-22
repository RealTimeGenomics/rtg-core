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
package com.rtg.vcf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import com.rtg.relation.Family;
import com.rtg.vcf.header.VcfHeader;

/**
 * Phases child genotype calls in VCF records according to pedigree. When a child genotype can be
 * unambiguously phased according to Mendelian inheritance, the genotype will be ordered such
 * that the allele inherited from the father is first, and the mothers is second.
 *
 */
public class ChildPhasingVcfAnnotator implements VcfAnnotator {

  private VcfHeader mHeader = null;

  private Collection<Family> mFamilies;

  /**
   * Constructor. All children within the supplied families will be phased with respect to their parents.
   * @param families the families of interest
   */
  public ChildPhasingVcfAnnotator(Family... families) {
    this(Arrays.asList(families));
  }

  /**
   * Constructor. All children within the supplied families will be phased with respect to their parents.
   * @param families the families of interest
   */
  public ChildPhasingVcfAnnotator(Collection<Family> families) {
    mFamilies = families;
  }

  @Override
  public void updateHeader(VcfHeader header) {
    mHeader = header;
  }

  @Override
  public void annotate(VcfRecord rec) {
    // Phase the calls
    final List<String> calls = rec.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE);
    final String[] phased = new String[calls.size()];
    for (Family f : mFamilies) { // The families should have at least two members with generated calls, but there may be other members missing
      final Integer fatherIndex = mHeader.getSampleIndex(f.getFather());
      final String fatherCall = (fatherIndex == null) ? VcfRecord.MISSING : calls.get(fatherIndex);
      final Integer motherIndex = mHeader.getSampleIndex(f.getMother());
      final String motherCall = (motherIndex == null) ? VcfRecord.MISSING : calls.get(motherIndex);
      for (String child : f.getChildren()) {
        final Integer childIndex = mHeader.getSampleIndex(child);
        if (childIndex != null) {
          phased[childIndex] = phaseDiploidCall(fatherCall, motherCall, calls.get(childIndex));
        }
      }
    }
    final ArrayList<String> newCalls = new ArrayList<>(calls.size());
    for (int i = 0; i < calls.size(); i++) {
      if (phased[i] == null) { // All not-yet phased calls get copied through unaltered
        phased[i] = calls.get(i);

        // If we wanted, we could phase homozygous remaining calls, but there is no information to be gained
      /*
      final int[] p = VcfUtils.splitGT(phased[i]);
      if (p.length == 2 && p[0] != -1 && p[0].equals(p[1])) {
        phased[i] = "" + p[0] + VcfUtils.PHASED_SEPARATOR + p[1];
      }
      */
      }
      newCalls.add(phased[i]);
    }
    // Plug phased calls into the record
    rec.getFormatAndSample().put(VcfUtils.FORMAT_GENOTYPE, newCalls);
  }

  /**
   * Attempt to phase the child with respect to the parents if it is non-missing and diploid
   * Assumes input obeys mendelian constraints.
   * If the call cannot be phased, returns null.
   * If the call can be phased, it is ordered such that the father allele comes first (this is for consistency with our child sim tools).
   */
  static String phaseDiploidCall(String fatherCall, String motherCall, String childCall) {
    final int[] childAlleles = VcfUtils.splitGt(childCall);
    if ((childAlleles.length == 2)
        && (childAlleles[0] != -1) && (childAlleles[1] != -1)) {

      if (childAlleles[0] == childAlleles[1]) {  // Homozygous always phased
        return "" + childAlleles[0] + VcfUtils.PHASED_SEPARATOR + childAlleles[1];
      } else {
        final int[] fatherAlleles = VcfUtils.splitGt(fatherCall);
        final int[] motherAlleles = VcfUtils.splitGt(motherCall);

        final boolean firstInFather = contains(fatherAlleles, childAlleles[0]);
        final boolean firstInMother = contains(motherAlleles, childAlleles[0]);
        final boolean secondInFather = contains(fatherAlleles, childAlleles[1]);
        final boolean secondInMother = contains(motherAlleles, childAlleles[1]);

        final boolean firstOnlyInFather = firstInFather && !firstInMother;
        final boolean firstOnlyInMother = firstInMother && !firstInFather;
        final boolean secondOnlyInFather = secondInFather && !secondInMother;
        final boolean secondOnlyInMother = secondInMother && !secondInFather;

        if ((firstOnlyInMother || secondOnlyInFather) && !(firstOnlyInFather || secondOnlyInMother)) {
          return "" + childAlleles[1] + VcfUtils.PHASED_SEPARATOR + childAlleles[0];
        } else if ((firstOnlyInFather || secondOnlyInMother) && !(firstOnlyInMother || secondOnlyInFather)) {
          return "" + childAlleles[0] + VcfUtils.PHASED_SEPARATOR + childAlleles[1];
        }
      }
    }
    return null;
  }

  private static boolean contains(int[] gt, int allele) {
    for (int current : gt) {
      if (allele == current) {
        return true;
      }
    }
    return false;
  }

}
