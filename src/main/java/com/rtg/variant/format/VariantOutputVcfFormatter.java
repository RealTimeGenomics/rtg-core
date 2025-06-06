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
package com.rtg.variant.format;

import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.EnumSet;

import com.rtg.relation.VcfPedigreeParser;
import com.rtg.util.MathUtils;
import com.rtg.util.PosteriorUtils;
import com.rtg.util.StringUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.util.VariantUtils;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

import htsjdk.samtools.SAMFileHeader;

/**
 * Handles the construction of an appropriate VcfHeader object, and conversion from a
 * Variant object into a VcfRecord object.
 *
 * Note: not thread safe.
 */
public class VariantOutputVcfFormatter {

  /** The default sample name */
  public static final String DEFAULT_SAMPLE = "SAMPLE";

  private final String[] mSampleNames;

  private final VariantParams mParams;

  private final EnumSet<VcfInfoField> mInfoFields = EnumSet.noneOf(VcfInfoField.class);
  private final EnumSet<VcfFilterField> mFilterFields = EnumSet.noneOf(VcfFilterField.class);
  private final EnumSet<VcfFormatField> mFormatFields = EnumSet.noneOf(VcfFormatField.class);

  /**
   * Create a VCF format variant output formatter
   * @param params variant output options params
   * @param sampleNames the sample names
   */
  public VariantOutputVcfFormatter(VariantParams params, String... sampleNames) {
    if (sampleNames.length == 0) {
      throw new IllegalArgumentException("No sample names");
    }
    if (sampleNames.length > 1 && params.genomeRelationships() == null) {
      throw new IllegalArgumentException("Too many sample names");
    }
    mSampleNames = Arrays.copyOf(sampleNames, sampleNames.length);
    mParams = params;
    initFields();
  }

  /**
   * Dummy constructor for tests or <code>toString</code> debug
   * @param sampleNames the sample names
   */
  public VariantOutputVcfFormatter(String... sampleNames) {
    mSampleNames = sampleNames.length > 0 ? sampleNames : new String[] {DEFAULT_SAMPLE};
    mParams = null;
    initFields();
  }

  private void initFields() {
    initInfoFields();
    initFilterFields();
    initFormatFields();
  }

  private void initInfoFields() {
    mInfoFields.add(VcfInfoField.XRX);
    mInfoFields.add(VcfInfoField.RCE);
    mInfoFields.add(VcfInfoField.CT);
    mInfoFields.add(VcfInfoField.DP_DPR);
    mInfoFields.add(VcfInfoField.NREF);
    if (mParams != null) {
      mInfoFields.addAll(mParams.infoAnnotations());
    }
  }

  private void initFilterFields() {
    mFilterFields.add(VcfFilterField.RC);
    mFilterFields.add(VcfFilterField.RX);
    mFilterFields.add(VcfFilterField.RCEQUIV);
    mFilterFields.add(VcfFilterField.PASS);
    if (mParams != null && mParams.regionsFilterBedFile() != null) {
      mFilterFields.add(VcfFilterField.BED);
    }
    mFilterFields.add(VcfFilterField.OTHER);
    mFilterFields.add(VcfFilterField.OC);
    if (mParams != null && mParams.maxAmbiguity() != null) {
      mFilterFields.add(VcfFilterField.A);
    }
    if (mParams != null) {
      if (mParams.ionTorrent()) {
        mFilterFields.add(VcfFilterField.IONT);
      }
    }
  }

  private void initFormatFields() {
    mFormatFields.add(VcfFormatField.GT);
    mFormatFields.add(VcfFormatField.DP);
    mFormatFields.add(VcfFormatField.AD);
    mFormatFields.add(VcfFormatField.GQ);
    mFormatFields.add(VcfFormatField.RE);
    mFormatFields.add(VcfFormatField.AR);
    mFormatFields.add(VcfFormatField.RS);
    mFormatFields.add(VcfFormatField.DPR);
    mFormatFields.add(VcfFormatField.ABP);
    mFormatFields.add(VcfFormatField.SBP);
    mFormatFields.add(VcfFormatField.RPB);
    mFormatFields.add(VcfFormatField.PPB);
    mFormatFields.add(VcfFormatField.PUR);
    mFormatFields.add(VcfFormatField.GL);
    mFormatFields.add(VcfFormatField.AQ);
    if (mParams != null) {
      if (mParams.vcfRp()) {
        mFormatFields.add(VcfFormatField.RP);
      }

      if (mParams.minVariantAllelicDepth() > 0 || mParams.minVariantAllelicFraction() > 0) {
        mFormatFields.add(VcfFormatField.VA);
        mFormatFields.add(VcfFormatField.ADE);
        mFormatFields.add(VcfFormatField.VAF);
      }
      mFormatFields.addAll(mParams.formatAnnotations());
    }
  }

  /**
   * Add these additional info fields to output.
   * @param extraInfoFields the additional info fields to add.
   */
  public void addExtraInfoFields(EnumSet<VcfInfoField> extraInfoFields) {
    mInfoFields.addAll(extraInfoFields);
  }

  /**
   * Add these additional info fields to output.
   * @param extraFormatFields the additional info fields to add.
   */
  public void addExtraFormatFields(EnumSet<VcfFormatField> extraFormatFields) {
    mFormatFields.addAll(extraFormatFields);
  }

  /**
   * Write VCF header to the output file. This is deprecated, you should instead use VcfWriter.
   * @param out where output is to be written.
   * @param params reified command line parameters.
   * @param header the SAM file header
   * @throws IOException when an IO error occurs
   */
  public void writeHeader(OutputStream out, VariantParams params, SAMFileHeader header) throws IOException {
    out.write(makeHeader(params, header).toString().getBytes());
  }

  /**
   * Format a variant call into VCF variant call output format. This is deprecated,
   * for writing to output files, instead use a VcfWriter.
   * @param call the call to format
   * @return the formatted string
   */
  public String formatCall(Variant call) {
    return makeVcfRecord(call).toString() + '\n';
  }


  /**
   * Create the VCF header corresponding to the current caller parameters.
   * @param params reified command line parameters.
   * @param header the SAM file header.
   * @return the VCF header.
   */
  public VcfHeader makeHeader(VariantParams params, SAMFileHeader header) {
    final VcfHeader vcf = new VcfHeader();
    vcf.addCommonHeader();
    //reference
    if (params.genome() != null) {
      vcf.addReference(params.genome().reader());
    }
    //contigs
    if (header != null) {
      vcf.addContigFields(header);
    }

    // Output SAMPLE / PEDIGREE lines
    if (params.genomeRelationships() != null) {
      VcfPedigreeParser.addPedigreeFields(vcf, params.genomeRelationships(), mSampleNames);
    }

    //work out sample header name
    for (final String sample : mSampleNames) {
      vcf.addSampleName(sample);
    }

    for (final VcfInfoField info : mInfoFields) {
      info.updateHeader(vcf);
    }

    for (final VcfFilterField filter : mFilterFields) {
      filter.updateHeader(vcf, mParams);
    }

    for (final VcfFormatField format : mFormatFields) {
      format.updateHeader(vcf);
    }

    return vcf;
  }

  /**
   * Format a variant call into VcfRecord call output format
   * @param call the call to format
   * @return the VcfRecord
   */
  public VcfRecord makeVcfRecord(Variant call) {
    //CHROM field

    //whether to include previous nt or not
    final boolean includePreviousNt = includePreviousNt(call);

    //REF field
    final StringBuilder finalRef = new StringBuilder();
    if (includePreviousNt) {
      final char previousRefNt = call.getLocus().getPreviousRefNt();
      finalRef.append(previousRefNt);
    }

    final String reference = call.getLocus().getRefNts();
    if (reference != null) {
      if (reference.length() > 0) {
        finalRef.append(reference);
      } else {
        assert includePreviousNt : "Not outputting reference, but also not outputting previous ref nt";
      }
    } else {
      finalRef.append('*');
    }

    final VcfRecord rec = new VcfRecord(call.getLocus().getSequenceName(), call.getLocus().getStart() - (includePreviousNt ? 1 : 0), finalRef.toString());

    //ID field
    rec.setId(VcfRecord.MISSING);

    rec.setNumberOfSamples(mSampleNames.length);
    //FILTER fields
    for (final VcfFilterField filter : mFilterFields) {
      filter.updateRecord(rec, call, mParams);
    }

    updateRecords(rec, call, includePreviousNt, false); // Add all annotations that operate off Variant/VariantSample
    updateRecords(rec, call, includePreviousNt, true);  // Add all annotations that operate off the partially constructed VcfRecord

    if (mSampleNames.length > 1) {
      setMulticallQuality(call, rec);
    } else {
      //ALT/QUAL field
      setCallQuality(call, rec);
    }

    return rec;
  }

  private void updateRecords(VcfRecord rec, Variant call, boolean includePreviousNt, boolean isVcfAnnotating) {
    //FORMAT fields
    //Determine fields which exist in at least one subcall and require the others to output a value also.
    final EnumSet<VcfFormatField> requiredFields = EnumSet.noneOf(VcfFormatField.class);
    for (final VcfFormatField format : mFormatFields) {
      if (format.isVcfAnnotator() == isVcfAnnotating) {
        for (int i = 0; i < mSampleNames.length; ++i) {
          if (format.hasValue(rec, call, call.getSample(i), mSampleNames[i], mParams)) {
            requiredFields.add(format);
            break;
          }
        }
      }
    }

    for (final VcfFormatField format : requiredFields) {
      format.updateRecord(rec, call, mSampleNames, mParams, includePreviousNt);
    }

    //INFO fields
    for (final VcfInfoField info : mInfoFields) {
      if (info.isVcfAnnotator() == isVcfAnnotating) {
        info.updateRecord(rec, call, mParams, includePreviousNt);
      }
    }
  }

  private static void setMulticallQuality(Variant call, VcfRecord rec) {
    if (call.getNonIdentityPosterior() == null) {
      rec.setQuality(VcfRecord.MISSING);
    } else if (rec.getAltCalls().isEmpty()) {
      rec.setQuality(MathUtils.cappedFloat(PosteriorUtils.nonIdentityPhredIfy(call.getNonIdentityPosterior())));
    } else {
      rec.setQuality(MathUtils.cappedFloat(PosteriorUtils.phredIfy(call.getNonIdentityPosterior())));
    }
  }

  private static void setCallQuality(Variant call, VcfRecord rec) {
    if (call.getSample(0).getName() != null) {
      if (call.getNonIdentityPosterior() == null) {
        rec.setQuality(VcfRecord.MISSING);
      } else if (rec.getAltCalls().isEmpty()) {
        rec.setQuality(MathUtils.cappedFloat(PosteriorUtils.nonIdentityPhredIfy(call.getNonIdentityPosterior())));
      } else {
        rec.setQuality(MathUtils.cappedFloat(PosteriorUtils.phredIfy(call.getNonIdentityPosterior())));
      }
    } else {
      rec.setQuality(VcfRecord.MISSING);
    }
  }

  /**
   * Check if the VCF record for a call would require the previous nucleotide to be included (usually for inserts/deletes)
   * @param call a variant call
   * @return true  if the record should prepend the previous nt
   */
  public static boolean includePreviousNt(Variant call) {
    boolean includePreviousNt = false;
    for (int i = 0; i < call.getNumberOfSamples(); ++i) {
      if (call.getSample(i) != null) {
        includePreviousNt |= includePreviousNt(call, call.getSample(i));
      }
    }
    return includePreviousNt;
  }

  private static boolean includePreviousNt(Variant v, VariantSample sample) {
    final String name = sample.getName();
    if (name == null) {
      return false; // this sample says nothing about whether to include
    }
    final String ref = v.getLocus().getRefNts();
    if (ref == null) {
      return false; // over coverage condition
    }
    if (ref.length() == 0) {
      return true;
    }
    if (checkCatLengths(name) || sample.getVariantAllele() != null && sample.getVariantAllele().length() == 0) {
      return true;
    }

    return false;
  }

  private static boolean checkCatLengths(String catName) {
    final String[] cats = StringUtils.split(catName, VariantUtils.COLON);
    for (final String c : cats) {
      if (c.length() == 0) {
        return true;
      }
    }
    return false;
  }
}
