/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
//  /** Mean normalized case coverage of bins in this segment. */
//  private static final String FORMAT_CASE = "NSC";
//  /** Mean normalized control coverage of bins in this segment. */
//  private static final String FORMAT_CTRL = "NCC";
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
    header.addInfoField(INFO_IMPRECISE, MetaType.FLAG, VcfNumber.FLAG, "Imprecise structural variation");
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
//    header.addFormatField(FORMAT_CASE, MetaType.FLOAT, VcfNumber.ONE, "Mean normalized case coverage of the segment");
//    header.addFormatField(FORMAT_CTRL, MetaType.FLOAT, VcfNumber.ONE, "Mean normalized control coverage of the segment");

    header.addSampleName(mSampleName);

    return header;
  }

  VcfRecord vcfRecord(Segment prev, Segment current, Segment next) throws IOException {
    // Spec says SV records are given the position BEFORE the SV.
    final int refPosition = Math.max(current.getStart() - 1, 0);
    final String ref = getRefBase(current.getSequenceName(), refPosition);
    final VcfRecord rec = new VcfRecord(current.getSequenceName(), refPosition, ref);

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
//    rec.addFormatAndSample(FORMAT_CASE, Utils.realFormat(current.meanNormalizedCaseCov(), 4));
//    rec.addFormatAndSample(FORMAT_CTRL, Utils.realFormat(current.meanNormalizedCtrlCov(), 4));

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
