/*
 * Copyright (c) 2017. Real Time Genomics Limited.
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

import static com.rtg.launcher.CommonFlags.FILE;
import static com.rtg.launcher.CommonFlags.INPUT_FLAG;
import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.rtg.alignment.SplitAlleles;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.SdfId;
import com.rtg.reader.SdfUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.Pair;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.vcf.eval.VcfEvalTask;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 * Attempt to decompose calls within VCF records into simpler forms.
 */
public final class VcfDecomposerCli extends AbstractCli {

  private static final String ORP = "ORP";
  private static final String ORL = "ORL";
  private long mCurrentSequenceId = -1;
  private byte[] mCurrentSequence = null;
  private long mTotalCallsSplit = 0;
  private long mTotalPieces = 0;
  private long mTotalRecords = 0;

  @Override
  public String moduleName() {
    return "vcfdecompose";
  }

  @Override
  public String description() {
    return "decompose complex variants within a VCF file";
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Decomposes complex variants within a VCF file.");
    CommonFlagCategories.setCategories(mFlags);
    mFlags.registerRequired('i', INPUT_FLAG, File.class, FILE, "VCF file containing variants to decompose. Use '-' to read from standard input").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT_FLAG, File.class, FILE, "output VCF file name. Use '-' to write to standard output").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('t', CommonFlags.TEMPLATE_FLAG, File.class, CommonFlags.SDF, "SDF of the reference genome the variants are called against").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(CommonFlags.NO_HEADER, "prevent VCF header from being written").setCategory(UTILITY);
    CommonFlags.initNoGzip(mFlags);
    CommonFlags.initIndexFlags(mFlags);
    CommonFlags.initForce(mFlags);
    mFlags.setValidator(new SnpAnnotatorValidator());
  }

  private static class SnpAnnotatorValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      return CommonFlags.validateInputFile(flags, INPUT_FLAG)
        && CommonFlags.validateOutputFile(flags, VcfUtils.getZippedVcfFileName(!flags.isSet(NO_GZIP), (File) flags.getValue(OUTPUT_FLAG)));
    }
  }

  private static void checkHeader(final VcfHeader header, final SdfId referenceSdfId) {
    final SdfId vcfSdfId = header.getSdfId();
    if (!vcfSdfId.check(referenceSdfId)) {
      Diagnostic.warning("Reference template ID mismatch, VCF variants were not created from the given reference");
    }
  }

  private void updateTemplate(final SequencesReader templateSequences, final long sequenceId) throws IOException {
    if (mCurrentSequenceId != sequenceId) {
      mCurrentSequenceId = sequenceId;
      mCurrentSequence = templateSequences.read(sequenceId);
    }
  }

  @Override
  protected int mainExec(final OutputStream out, final PrintStream err) throws IOException {
    final File inputFile = (File) mFlags.getValue(INPUT_FLAG);
    final File output = (File) mFlags.getValue(OUTPUT_FLAG);
    final boolean gzip = !mFlags.isSet(NO_GZIP);
    final boolean stdout = FileUtils.isStdio(output);
    final File templateFile = (File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG);
    SdfUtils.validateHasNames(templateFile);
    try (final SequencesReader templateSequences = SequencesReaderFactory.createDefaultSequencesReader(templateFile, LongRange.NONE)) {
      SdfUtils.validateNoDuplicates(templateSequences, false);
      final Map<String, Long> nameMap = VcfEvalTask.getNameMap(templateSequences);
      try (VcfReader reader = VcfReader.openVcfReader(inputFile)) {
        final VcfHeader header = reader.getHeader();
        checkHeader(header, templateSequences.getSdfId());
        header.ensureContains(new InfoField(ORP, MetaType.STRING, VcfNumber.ONE, "Original position"));
        header.ensureContains(new InfoField(ORL, MetaType.STRING, VcfNumber.ONE, "Original longest allele"));
        final File vcfFile = stdout ? null : VcfUtils.getZippedVcfFileName(gzip, output);
        try (VcfWriter writer = new VcfWriterFactory(mFlags).addRunInfo(true).make(header, vcfFile, out)) {
          while (reader.hasNext()) {
            final VcfRecord rec = reader.next();
            updateTemplate(templateSequences, nameMap.get(rec.getSequenceName()));
            for (final VcfRecord res : decompose(rec)) {
              writer.write(res);
            }
          }
        }
      }
    }
    try (final PrintStream ps = new PrintStream(out)) {
      ps.println("Total records : " + mTotalRecords);
      ps.println("Number of records decomposed : " + mTotalCallsSplit);
      ps.println("Remaining records : " + (mTotalRecords + (mTotalPieces - mTotalCallsSplit)));
    }
    return 0;
  }

  private boolean needsAnchorBase(final String[] alleles) {
    for (final String a : alleles) {
      if (a.isEmpty()) {
        return true;
      }
    }
    return false;
  }

  private int longestAllele(final String ref, final List<String> alts) {
    int length = ref.length();
    for (final String a : alts) {
      length = Math.max(length, a.length());
    }
    return length;
  }

  private String getAnchoredAllele(final String allele, final int start, final boolean needAnchor) {
    return needAnchor ? DnaUtils.base(mCurrentSequence, start - 1) + allele : allele;
  }

  private int[] updateAlts(final List<String> altCalls, final String[] alleles, final int start, final boolean needAnchor) {
    final int[] altRemap = new int[alleles.length]; // Map old allele number to new allele number
    final String refAllele = getAnchoredAllele(alleles[0], start, needAnchor);
    altCalls.clear();
    for (int k = 1; k < alleles.length; ++k) {
      final String allele = alleles[k];
      if (allele.equals(refAllele)) {
        altRemap[k] = 0; // This allele is now reference
      } else {
        final String newAllele = getAnchoredAllele(allele, start, needAnchor);
        final int existing = altCalls.indexOf(newAllele);
        if (existing >= 0) {
          // We already have this alternative
          altRemap[k] = existing + 1;
        } else {
          // This is a new alternative
          altCalls.add(newAllele);
          altRemap[k] = altCalls.size();
        }
      }
    }
    return altRemap;
  }

  private static String makeGt(final int[] gt, boolean isPhased) {
    final StringBuilder sb = new StringBuilder();
    for (final int a : gt) {
      if (sb.length() > 0) {
        sb.append(isPhased ? VcfUtils.PHASED_SEPARATOR : VcfUtils.UNPHASED_SEPARATOR);
      }
      if (a == VcfUtils.MISSING_GT) {
        sb.append(VcfUtils.MISSING_VALUE);
      } else {
        sb.append(a);
      }
    }
    return sb.toString();
  }

  private List<VcfRecord> decompose(final VcfRecord rec) {
    ++mTotalRecords;
    final SplitAlleles pwa = SplitAlleles.create(rec);
    final List<Pair<Integer, String[]>> split = SplitAlleles.removeAllRef(pwa.partition());
    if (split.size() <= 1) {
      return Collections.singletonList(rec); // Efficiency, no change in record
    }
    ++mTotalCallsSplit;
    final ArrayList<VcfRecord> res = new ArrayList<>(split.size());
    for (final Pair<Integer, String[]> s : split) {
      final VcfRecord splitRecord = new VcfRecord(rec);
      splitRecord.addInfo("ORP", String.valueOf(rec.getStart() + 1)); // 1-based for output
      splitRecord.addInfo("ORL", String.valueOf(longestAllele(rec.getRefCall(), rec.getAltCalls())));
      final String[] alleles = s.getB();
      final boolean needAnchor = needsAnchorBase(alleles);
      final Integer offset = s.getA();
      final int start = splitRecord.getStart() + offset;
      if (needAnchor) {
        splitRecord.setStart(start - 1);
        splitRecord.setRefCall(getAnchoredAllele(alleles[0], start, true));
      } else {
        splitRecord.setStart(start);
        splitRecord.setRefCall(alleles[0]);
      }
      final int[] alleleMap = updateAlts(splitRecord.getAltCalls(), alleles, start, needAnchor);
      // Go through samples and update genotypes
      final List<String> oldGenotypes = splitRecord.getFormat(VcfUtils.FORMAT_GENOTYPE);
      for (int sample = 0; sample < splitRecord.getNumberOfSamples(); ++sample) {
        final String oldGt = oldGenotypes.get(sample);
        final int[] gt = VcfUtils.splitGt(oldGt);
        for (int k = 0; k < gt.length; ++k) {
          gt[k] = gt[k] == VcfUtils.MISSING_GT ? VcfUtils.MISSING_GT : alleleMap[gt[k]];
        }
        final boolean isPhased = oldGt.indexOf(VcfUtils.PHASED_SEPARATOR) >= 0;
        splitRecord.setFormatAndSample(VcfUtils.FORMAT_GENOTYPE, makeGt(gt, isPhased), sample);
      }
      ++mTotalPieces;
      res.add(splitRecord);
    }
    return res;
  }

}
