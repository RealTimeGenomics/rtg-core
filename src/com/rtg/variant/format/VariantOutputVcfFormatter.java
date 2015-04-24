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
package com.rtg.variant.format;

import java.io.IOException;
import java.io.OutputStream;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;

import com.rtg.relation.VcfPedigreeParser;
import com.rtg.util.MathUtils;
import com.rtg.util.StringUtils;
import com.rtg.variant.PosteriorUtils;
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
 * Note: not multi thread safe.
 */
public class VariantOutputVcfFormatter {

  /** The default sample name */
  public static final String DEFAULT_SAMPLE = "SAMPLE";

  private final Map<String, Integer> mSampleColumns;
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
    mSampleNames = new String[sampleNames.length];
    if (sampleNames.length > 1 && params.genomeRelationships() == null) {
      throw new IllegalArgumentException("Too many sample names");
    }
    mSampleColumns = new HashMap<>();
    for (int i = 0; i < mSampleNames.length; i++) {
      mSampleNames[i] = sampleNames[i];
      mSampleColumns.put(mSampleNames[i], i);
    }
    mParams = params;
    initFields();
  }

  /**
   * Dummy constructor for tests or <code>toString</code> debug
   * @param sampleNames the sample names
   */
  public VariantOutputVcfFormatter(String... sampleNames) {
    mSampleNames = sampleNames.length > 0 ? sampleNames : new String[] {DEFAULT_SAMPLE};
    mSampleColumns = null;
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
    mInfoFields.add(VcfInfoField.DPR);
    if (mSampleNames.length > 1) {
      mInfoFields.add(VcfInfoField.DP);
    }
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
    if (mSampleNames.length > 1 && mSampleColumns != null) {
      mFormatFields.add(VcfFormatField.RQ);
      mFormatFields.add(VcfFormatField.DN);
      mFormatFields.add(VcfFormatField.DNP);
    }
    if (mParams != null) {
      if (mParams.vcfRp()) {
        mFormatFields.add(VcfFormatField.RP);
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
    final VcfRecord rec = new VcfRecord();

    rec.setNumberOfSamples(mSampleNames.length);
    //CHROM field
    rec.setSequence(call.getLocus().getSequenceName());

    //whether to include previous nt or not
    final boolean includePreviousNt = includePreviousNt(call);

    rec.setStart(call.getLocus().getStart() - (includePreviousNt ? 1 : 0));

    //ID field
    rec.setId(VcfRecord.MISSING);

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
    rec.setRefCall(finalRef.toString());

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
        for (int i = 0; i < mSampleNames.length; i++) {
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
    } else if (rec.getAltCalls().size() == 0) {
      rec.setQuality(MathUtils.cappedFloat(PosteriorUtils.nonIdentityPhredIfy(call.getNonIdentityPosterior())));
    } else {
      rec.setQuality(MathUtils.cappedFloat(PosteriorUtils.phredIfy(call.getNonIdentityPosterior())));
    }
  }

  private static void setCallQuality(Variant call, VcfRecord rec) {
    if (call.getSample(0).getName() != null) {
      if (call.getNonIdentityPosterior() == null) {
        rec.setQuality(VcfRecord.MISSING);
      } else if (rec.getAltCalls().size() == 0) {
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
    for (int i = 0; i < call.getNumberOfSamples(); i++) {
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
    if (checkCatLengths(name)) {
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
