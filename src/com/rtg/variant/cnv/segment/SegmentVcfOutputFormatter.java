/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.cnv.segment;

import static com.rtg.vcf.VcfUtils.ALT_DEL;
import static com.rtg.vcf.VcfUtils.ALT_DUP;
import static com.rtg.vcf.VcfUtils.FORMAT_GENOTYPE;
import static com.rtg.vcf.VcfUtils.INFO_CIEND;
import static com.rtg.vcf.VcfUtils.INFO_CIPOS;
import static com.rtg.vcf.VcfUtils.INFO_END;
import static com.rtg.vcf.VcfUtils.INFO_IMPRECISE;
import static com.rtg.vcf.VcfUtils.INFO_SVTYPE;
import static com.rtg.vcf.VcfUtils.SvType.DEL;
import static com.rtg.vcf.VcfUtils.SvType.DUP;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import com.rtg.mode.DnaUtils;
import com.rtg.reader.NamesInterface;
import com.rtg.reader.SequencesReader;
import com.rtg.util.Utils;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.AltField;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 * Formats segment output as VCF
 */
public class SegmentVcfOutputFormatter {

  // Bin count
  private static final String INFO_BC = "BC";
  // Read depth ratio
  static final String FORMAT_RDR = "RDR";
  // Log of ratio
  static final String FORMAT_LOGR = "LR";
  /** The FORMAT field we use to store an overall quality score */
  public static final String FORMAT_SQS = "SQS";


  private final SequencesReader mTemplate;
  private final double mThreshold;
  private final Map<String, Integer> mSequenceMap = new HashMap<>();
  private final int mMinBins;
  private final String mSampleName;


  private String mCurrentSequenceName;
  private int mCurrentSequenceId;
  private int mCurrentSequenceLength;

  /**
   * Create a new object
   * @param genomeSequences reader for template
   * @param threshold applied to log ratio beyond which a copy number alteration is output
   * @param minBins minimum number of bins before a copy number alteration is called
   * @param sampleName the name to use for the sample column
   * @throws IOException if error
   */
  public SegmentVcfOutputFormatter(SequencesReader genomeSequences, double threshold, int minBins, String sampleName) throws IOException {
    mTemplate = genomeSequences;
    mThreshold = threshold;
    mMinBins = minBins;
    mSampleName = sampleName;
    final NamesInterface pni = genomeSequences.names();
    for (long i = 0; i < pni.length(); ++i) {
      mSequenceMap.put(genomeSequences.names().name(i), (int) i);
    }
    assert mThreshold >= 0;
  }

  VcfHeader header() throws IOException {
    final VcfHeader header = new VcfHeader();
    header.addCommonHeader();
    header.addReference(mTemplate);
    header.addContigFields(mTemplate);
    header.addAltField(new AltField(DEL.name(), "Deletion"));
    header.addAltField(new AltField(DUP.name(), "Duplication"));

    header.addInfoField(INFO_END, MetaType.INTEGER, VcfNumber.ONE, "End position of the variant described in this record");
    header.addInfoField(INFO_IMPRECISE, MetaType.FLAG, new VcfNumber("0"), "Imprecise structural variation");
    header.addInfoField(INFO_SVTYPE, MetaType.STRING, VcfNumber.ONE, "Type of structural variant");
    header.addInfoField(INFO_BC, MetaType.INTEGER, VcfNumber.ONE, "Number of bins contained within the region");
    header.addInfoField(INFO_CIPOS, MetaType.INTEGER, new VcfNumber("2"), "Confidence interval around POS for imprecise variants");
    header.addInfoField(INFO_CIEND, MetaType.INTEGER, new VcfNumber("2"), "Confidence interval around END for imprecise variants");
    //Region/gene names contained within each segment?
    //header.addInfoField(INFO_GENE, MetaType.STRING, VcfNumber.ONE, "Bin names");

    VcfFormatField.GT.updateHeader(header);
    header.addFormatField(FORMAT_SQS, MetaType.FLOAT, VcfNumber.ONE, "Segment quality score");
    header.addFormatField(FORMAT_RDR, MetaType.FLOAT, VcfNumber.ONE, "Mean normalized RD ratio with respect to control");
    header.addFormatField(FORMAT_LOGR, MetaType.FLOAT, VcfNumber.ONE, "Log2 of RD ratio with respect to control");

    header.addSampleName(mSampleName);

    return header;
  }

  VcfRecord vcfRecord(String seqName, Segment prev, Segment current, Segment next) throws IOException {
    // Spec says SV records are given the position BEFORE the SV.
    final int refPosition = Math.max(current.getStart() - 1, 0);
    final String ref = getRefBase(seqName, refPosition);
    final VcfRecord rec = new VcfRecord(seqName, refPosition, ref);

    // Classify the segment as gain or loss
    boolean altered = false;
    if (current.bins() >= mMinBins) {
      if (current.mean() >= mThreshold) {
        altered = true;
        rec.addAltCall(ALT_DUP);
        rec.setInfo(INFO_SVTYPE, DUP.name());
      } else if (current.mean() <= -mThreshold) {
        altered = true;
        rec.addAltCall(ALT_DEL);
        rec.setInfo(INFO_SVTYPE, DEL.name());
      }
    }

    rec.setInfo(INFO_END, String.valueOf(current.getEnd()));
    rec.setInfo(INFO_IMPRECISE);
    // For CIPOS and CIEND, we extend outward to the closest neighboring segment boundary (or the edge of the reference sequence)
    // For the inward side, we go to half of the bin length
    // Thus the extended CI span of one segment will overlap that of the next.
    final int pend = prev == null ? 0 : prev.getEnd();
    final String cipos = (pend - current.getStart())
      + ","
      + (current.firstBinLength() / 2);
    rec.setInfo(INFO_CIPOS, cipos);

    final int nstart = next == null ? mCurrentSequenceLength : next.getStart();
    final String ciend = (-current.lastBinLength() / 2)
      + ","
      + (nstart - current.getEnd());
    rec.setInfo(INFO_CIEND, ciend);

    rec.setInfo(INFO_BC, String.valueOf(current.bins()));

    rec.setNumberOfSamples(1);
    rec.addFormatAndSample(FORMAT_GENOTYPE, altered ? "1" : "0");
    rec.addFormatAndSample(FORMAT_LOGR, Utils.realFormat(current.mean(), 4));
    rec.addFormatAndSample(FORMAT_RDR, Utils.realFormat(rdr(current.mean()), 4));
    rec.addFormatAndSample(FORMAT_SQS, Utils.realFormat(sqs(current.mean()), 4)); // For now, just use abs of LogR as proxy for quality

    return rec;
  }

  static double rdr(final double lr) {
    return Math.pow(2, lr);
  }

  static double sqs(final double lr) {
    return Math.abs(lr);
  }

  private String getRefBase(String refName, int pos) throws IOException {
    if (!refName.equals(mCurrentSequenceName)) {
      final Integer id = mSequenceMap.get(refName);
      if (id == null) {
        throw new RuntimeException("Reference sequence '" + refName + "' was not contained in supplied SDF");
      }
      mCurrentSequenceName = refName;
      mCurrentSequenceId = id;
      mCurrentSequenceLength = mTemplate.length(id);
    }
    assert pos < mCurrentSequenceLength : "pos=" + pos + " currentLen=" + mCurrentSequenceLength;
    if (pos < 0 || (pos > mCurrentSequenceLength - 1)) {
      return "N";
    }
    final byte[] temp = new byte[1];
    mTemplate.read(mCurrentSequenceId, temp, pos, 1);
    return DnaUtils.bytesToSequenceIncCG(temp);
  }

}
