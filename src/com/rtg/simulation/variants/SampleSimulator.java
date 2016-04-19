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
package com.rtg.simulation.variants;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.reader.SequencesReader;
import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.ReferenceGenome.ReferencePloidy;
import com.rtg.reference.ReferenceSequence;
import com.rtg.reference.Sex;
import com.rtg.simulation.SimulationUtils;
import com.rtg.util.PortableRandom;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.io.FileUtils;
import com.rtg.vcf.VariantStatistics;
import com.rtg.vcf.VariantType;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.VcfWriter;
import com.rtg.vcf.header.FormatField;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 * Creates a genotyped sample using population variants defined in a VCF file.
 */
public class SampleSimulator {

  protected final SequencesReader mReference;
  private final PortableRandom mRandom;
  private final ReferencePloidy mDefaultPloidy;
  private VariantStatistics mStats = null;
  private boolean mSeenVariants = false;
  private int mDefaultAfCount = 0;

  /**
   * @param reference input reference data
   * @param rand random number generator
   * @param ploidy the default ploidy to use if no reference specification is present
   */
  public SampleSimulator(SequencesReader reference, PortableRandom rand, ReferencePloidy ploidy) {
    mReference = reference;
    mRandom = rand;
    mDefaultPloidy = ploidy;
  }

  private static VcfHeader getVcfHeader(File vcfFile) throws IOException {
    try (VcfReader reader = VcfReader.openVcfReader(vcfFile)) {
      return reader.getHeader();
    }
  }

  /**
   * Create a genotyped sample using population variants defined in file.
   * @param vcfPopFile input population data. requires allele frequencies
   * @param vcfOutFile destination of sample genotype
   * @param sample name to give the generated sample
   * @param sex sex of the generated sample
   * @throws java.io.IOException if an IO error occurs
   */
  public void mutateIndividual(File vcfPopFile, File vcfOutFile, String sample, Sex sex) throws IOException {
    final VcfHeader header = getVcfHeader(vcfPopFile);
    if (header.getSampleNames().contains(sample)) {
      throw new NoTalkbackSlimException("sample '" + sample + "' already exists");
    }
    header.addSampleName(sample);
    mStats = new VariantStatistics(null);
    mStats.onlySamples(sample);
    boolean foundGt = false;
    for (FormatField ff : header.getFormatLines()) {
      if (ff.getId().equals(VcfUtils.FORMAT_GENOTYPE)) {
        foundGt = true;
        break;
      }
    }
    if (!foundGt) {
      header.addFormatField(VcfUtils.FORMAT_GENOTYPE, MetaType.STRING, VcfNumber.ONE, "Genotype");
    }
    if (sex == Sex.FEMALE || sex == Sex.MALE) {
      header.addLine(VcfHeader.SAMPLE_STRING + "=<ID=" + sample + ",Sex=" + sex.toString() + ">");
    }
    header.addRunInfo();
    header.addLine(VcfHeader.META_STRING + "SEED=" + mRandom.getSeed());

    try (VcfWriter vcfOut = new VcfWriter(header, vcfOutFile, null, FileUtils.isGzipFilename(vcfOutFile), true)) {
      final ReferenceGenome refG = new ReferenceGenome(mReference, sex, mDefaultPloidy);
      for (long i = 0; i < mReference.numberSequences(); i++) {
        final ReferenceSequence refSeq = refG.sequence(mReference.name(i));
        mutateSequence(vcfPopFile, vcfOut, refSeq);
      }
    }
    if (!mSeenVariants) {
      Diagnostic.warning("No input variants! (is the VCF empty, or against an incorrect reference?)");
      Diagnostic.info("");
    } else if (mDefaultAfCount > 0) {
      Diagnostic.warning("" + mDefaultAfCount + " input records had no allele frequency information.");
      Diagnostic.info("");
    }

    Diagnostic.info(mStats.getStatistics());
  }

  //writes sample to given writer, returns records as list
  private List<VcfRecord> mutateSequence(File vcfPopFile, VcfWriter vcfOut, ReferenceSequence refSeq) throws IOException {
    Diagnostic.userLog("Simulating mutations on sequence: " + refSeq.name());
    final ArrayList<VcfRecord> sequenceMutations = new ArrayList<>();
    final int ploidyCount = refSeq.ploidy().count() >= 0 ? refSeq.ploidy().count() : 1; //effectively treats polyploid as haploid
    try (VcfReader reader = VcfReader.openVcfReader(vcfPopFile, new RegionRestriction(refSeq.name()))) {
      int lastVariantEnd = -1;
      while (reader.hasNext()) {
        mSeenVariants = true;
        final VcfRecord v = reader.next();
        v.addFormat(VcfUtils.FORMAT_GENOTYPE); // Ensure the record has a notion of genotype - strictly we should also initialize any pre-exising sample columns with empty GT
        final StringBuilder gt = new StringBuilder();
        if (refSeq.ploidy().count() != 0) {
          final List<String> allFreqStr = v.getInfo().get(VcfUtils.INFO_ALLELE_FREQ);
          final double[] dist;
          if (allFreqStr == null) {
            mDefaultAfCount++;
            final double[] defaultDist = new double[v.getAltCalls().size() + 1];
            Arrays.fill(defaultDist, 1.0); // All alleles are equally likely
            dist = SimulationUtils.cumulativeDistribution(defaultDist);
          } else {
            dist = new double[allFreqStr.size() + 1];
            double ac = 0.0;
            for (int i = 0; i < allFreqStr.size(); i++) {
              ac += Double.parseDouble(allFreqStr.get(i));
              dist[i] = ac;
            }
            if (ac > 1.0) {
              throw new NoTalkbackSlimException("Sum of AF probabilities exceeds 1.0 for record " + v);
            }
          }
          dist[dist.length - 1] = 1.0; //the reference consumes the rest of the distribution

          int variantEnd = lastVariantEnd; // Need to use separate variable while going through the ploidy loop, otherwise we won't make hom variants.
          for (int i = 0; i < ploidyCount; i++) {
            if (gt.length() != 0) {
              gt.append(VcfUtils.PHASED_SEPARATOR);
            }
            final int alleleId;
            if (v.getStart() < lastVariantEnd) { // Do not generate overlapping variants
              alleleId = 0;
            } else {
              final double d = mRandom.nextDouble();
              alleleId = chooseAllele(dist, d, v.getAltCalls().size());
            }
            final VariantType svType = alleleId == 0 ? null : VariantType.getSymbolicAlleleType(v.getAltCalls().get(alleleId - 1));
            if (svType != null) {
              Diagnostic.warning("Symbolic variant ignored: " + v);
              variantEnd = Math.max(variantEnd, v.getEnd());
              gt.append(0); // For now, skip symbolic alleles.
            } else {
              // Updating variantEnd even for alleleId == 0 prevents VCF representational issues with overlapping variants.
              //if (alleleId > 0) {
              variantEnd = Math.max(variantEnd, v.getEnd());
              //}
              gt.append(alleleId);
            }
          }
          lastVariantEnd = variantEnd;
        } else {
          //ploidy count == 0
          gt.append('.');
        }
        v.setNumberOfSamples(v.getNumberOfSamples() + 1);
        for (String format : v.getFormats()) {
          final String value = format.equals(VcfUtils.FORMAT_GENOTYPE) ? gt.toString() : VcfRecord.MISSING;
          v.addFormatAndSample(format, value);
        }
        sequenceMutations.add(v);
        vcfOut.write(v);
        mStats.tallyVariant(vcfOut.getHeader(), v);
      }
    }
    return sequenceMutations;
  }

  private static int chooseAllele(double[] dist, double d, int size) {
    for (int j = 0; j < dist.length; j++) {
      if (d < dist[j]) {
        if (j < size) {
          return j + 1;
        } else {
          return 0;
        }
      }
    }
    throw new RuntimeException("d: " + d + " should be less than dist[length - 1]:" +  dist[dist.length - 1]);
  }
}
