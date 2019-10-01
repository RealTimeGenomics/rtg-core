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

package com.rtg.variant.avr;

import static com.rtg.launcher.CommonFlags.BED_REGIONS_FLAG;
import static com.rtg.launcher.CommonFlags.FILE;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Properties;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.sam.SamRangeUtils;
import com.rtg.util.ThreadAware;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.annotation.DerivedAnnotations;
import com.rtg.vcf.header.VcfHeader;

/**
 */
public class BuilderCli extends AbstractCli {

  private static final String OUTPUT_FLAG = "output";
  private static final String POSITIVE_VCF_FLAG = "positive";
  private static final String NEGATIVE_VCF_FLAG = "negative";
  private static final String ANNOTATED_VCF_FLAG = "annotated";
  private static final String SAMPLE_FLAG = "sample";
  private static final String INFO_ANNOTATIONS_FLAG = "info-annotations";
  private static final String FORMAT_ANNOTATIONS_FLAG = "format-annotations";
  private static final String DERIVED_ANNOTATIONS_FLAG = "derived-annotations";
  private static final String QUAL_ANNOTATION_FLAG = "qual-annotation";
  private static final String X_MODEL_PARAMS_FLAG = "Xmodel-params";
  private static final String X_MODEL_TYPE_FLAG = "Xmodel-type";
  private static final String X_DUMP_MODEL_FLAG = "Xdump-model";
  private static final String X_POS_WEIGHT = "Xpos-weight";
  private static final String X_NEG_WEIGHT = "Xneg-weight";

  @Override
  public String moduleName() {
    return "avrbuild";
  }

  @Override
  public String description() {
    return "build an AVR model from training examples";
  }

  @Override
  protected void initFlags() {
    CommonFlagCategories.setCategories(mFlags);
    mFlags.setDescription("Create an AVR model from positive and negative training examples.");

    mFlags.registerOptional(BED_REGIONS_FLAG, File.class, FILE, "if set, only read VCF records that overlap the ranges contained in the specified BED file").setCategory(INPUT_OUTPUT);

    mFlags.registerRequired('o', OUTPUT_FLAG, File.class, FILE, "output AVR model").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional('p', POSITIVE_VCF_FLAG, File.class, FILE, "VCF file containing positive training examples").setCategory(INPUT_OUTPUT).setMaxCount(Integer.MAX_VALUE);
    mFlags.registerOptional('n', NEGATIVE_VCF_FLAG, File.class, FILE, "VCF file containing negative training examples").setCategory(INPUT_OUTPUT).setMaxCount(Integer.MAX_VALUE);
    mFlags.registerOptional('a', ANNOTATED_VCF_FLAG, File.class, FILE, "VCF file containing training examples annotated with CALL=TP/FP").setCategory(INPUT_OUTPUT).setMaxCount(Integer.MAX_VALUE);

    mFlags.registerOptional('s', SAMPLE_FLAG, String.class, CommonFlags.STRING, "the name of the sample to select (required when using multi-sample VCF files)").setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(INFO_ANNOTATIONS_FLAG, String.class, CommonFlags.STRING, "INFO fields to use in model").setCategory(SENSITIVITY_TUNING).setMaxCount(Integer.MAX_VALUE).enableCsv();
    mFlags.registerOptional(FORMAT_ANNOTATIONS_FLAG, String.class, CommonFlags.STRING, "FORMAT fields to use in model").setCategory(SENSITIVITY_TUNING).setMaxCount(Integer.MAX_VALUE).enableCsv();
    mFlags.registerOptional(QUAL_ANNOTATION_FLAG, "if set, use QUAL annotation in model").setCategory(SENSITIVITY_TUNING);
    final EnumSet<DerivedAnnotations> derivedAnnotations = DerivedAnnotations.singleValueAnnotations();
    final List<String> derivedRange = new ArrayList<>(derivedAnnotations.size());
    for (final DerivedAnnotations derived : derivedAnnotations) {
      derivedRange.add(derived.toString());
    }
    mFlags.registerOptional(DERIVED_ANNOTATIONS_FLAG, String.class, CommonFlags.STRING, "derived fields to use in model").setParameterRange(derivedRange).setMaxCount(Integer.MAX_VALUE).enableCsv().setCategory(SENSITIVITY_TUNING);

    mFlags.registerOptional(X_MODEL_PARAMS_FLAG, File.class, "PROPERTIES", "property file containing model parameters").setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(X_POS_WEIGHT, Double.class, CommonFlags.FLOAT, "weight to assign to positive training examples", 1.0).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(X_NEG_WEIGHT, Double.class, CommonFlags.FLOAT, "weight to assign to negative training examples", 1.0).setCategory(SENSITIVITY_TUNING);

    CommonFlags.initThreadsFlag(mFlags);
    mFlags.registerOptional(X_MODEL_TYPE_FLAG, ModelType.class, "TYPE", "specify the type of AVR model to build", ModelType.ML).setCategory(UTILITY);
    mFlags.registerOptional(X_DUMP_MODEL_FLAG, "dump model contents").setCategory(INPUT_OUTPUT);
    mFlags.setValidator(flags -> flags.checkOr(INFO_ANNOTATIONS_FLAG, FORMAT_ANNOTATIONS_FLAG, QUAL_ANNOTATION_FLAG, DERIVED_ANNOTATIONS_FLAG)
      && checkTrainingLabels(flags));
  }

  private boolean checkTrainingLabels(CFlags flags) {
    if (!flags.isSet(ANNOTATED_VCF_FLAG) && !(flags.isSet(POSITIVE_VCF_FLAG) && flags.isSet(NEGATIVE_VCF_FLAG))) {
      flags.setParseMessage("When not using annotated training VCFs, both positive and negative training VCFs are required");
      return false;
    }
    return true;
  }


  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {

    final long startMillis = System.currentTimeMillis();
    final File outFile = (File) mFlags.getValue(OUTPUT_FLAG);
    final String[] infoAttributes = getStrings(INFO_ANNOTATIONS_FLAG);
    final String[] formatAttributes = getStrings(FORMAT_ANNOTATIONS_FLAG);
    final String[] derivedAttributes = getStrings(DERIVED_ANNOTATIONS_FLAG);

    final AbstractModelBuilder<?> builder;
    final ModelType modelType = (ModelType) mFlags.getValue(X_MODEL_TYPE_FLAG);
    switch (modelType) {
      case ML:
        builder = new MlAvrModelBuilder(formatAttributes, infoAttributes, derivedAttributes);
        break;
      case GT_COMPLEX:
        builder = new GtQualComplexMultiplierModelBuilder(formatAttributes, infoAttributes, derivedAttributes);
        break;
      case ARFF:
        builder = new ArffDatasetWriter(formatAttributes, infoAttributes, derivedAttributes);
        break;
      case NULL:
        builder = new NullModelBuilder(formatAttributes, infoAttributes, derivedAttributes);
        break;
      default:
        throw new IllegalArgumentException("Unhandled model type: " + modelType);
    }
    builder.useQualAttribute(mFlags.isSet(QUAL_ANNOTATION_FLAG));
    if (builder instanceof ThreadAware) {
      ((ThreadAware) builder).setNumberOfThreads(CommonFlags.parseThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG)));
    }
    if (mFlags.isSet(X_MODEL_PARAMS_FLAG)) {
      final Properties parameters = new Properties();
      try (FileInputStream fis = new FileInputStream((File) mFlags.getValue(X_MODEL_PARAMS_FLAG))) {
        parameters.load(fis);
      }
      builder.setModelParameters(parameters);
    }

    final ReferenceRanges<String> ranges;
    if (mFlags.isSet(BED_REGIONS_FLAG)) {
      Diagnostic.developerLog("Loading BED regions");
      ranges = SamRangeUtils.createBedReferenceRanges((File) mFlags.getValue(BED_REGIONS_FLAG));
    } else {
      ranges = null;
    }

    final String sampleName = (String) mFlags.getValue(SAMPLE_FLAG);

    final ArrayList<VcfDataset> datasets = new ArrayList<>();
    for (Object o : mFlags.getValues(POSITIVE_VCF_FLAG)) {
      final File posVcf = (File) o;
      datasets.add(new VcfDataset(posVcf, getSampleNumber(posVcf, sampleName), VcfDataset.Classifications.ALL_POSITIVE, !mFlags.isSet(X_POS_WEIGHT), (Double) mFlags.getValue(X_POS_WEIGHT), ranges));
    }
    for (Object o : mFlags.getValues(NEGATIVE_VCF_FLAG)) {
      final File negVcf = (File) o;
      datasets.add(new VcfDataset(negVcf, getSampleNumber(negVcf, sampleName), VcfDataset.Classifications.ALL_NEGATIVE,  !mFlags.isSet(X_NEG_WEIGHT), (Double) mFlags.getValue(X_NEG_WEIGHT), ranges));
    }
    for (Object o : mFlags.getValues(ANNOTATED_VCF_FLAG)) {
      final File aVcf = (File) o;
      datasets.add(new VcfDataset(aVcf, getSampleNumber(aVcf, sampleName), VcfDataset.Classifications.ANNOTATED, !mFlags.isSet(X_NEG_WEIGHT), (Double) mFlags.getValue(X_NEG_WEIGHT), ranges));
    }

    try {
      builder.build(datasets.toArray(new VcfDataset[0]));
      builder.save(outFile);
    } catch (IllegalArgumentException iae) {
      throw new NoTalkbackSlimException(iae.getMessage());
    }

    try (final PrintStream ps = new PrintStream(out)) {
      if (mFlags.isSet(X_DUMP_MODEL_FLAG)) {
        ps.println(builder.getModel().toString());
      }
      ps.println(LoggedCli.getDuration(startMillis, true));
    }
    return 0;
  }

  private String[] getStrings(String flag) {
    final String[] derivedAttributes;
    if (mFlags.isSet(flag)) {
      derivedAttributes = mFlags.getValues(flag).stream().map(Object::toString).toArray(String[]::new);
    } else {
      derivedAttributes = new String[0];
    }
    return derivedAttributes;
  }

  private static int getSampleNumber(File vcf, String sampleName) throws IOException {
    final VcfHeader header;
    try (final VcfReader reader = VcfReader.openVcfReader(vcf)) {
      header = reader.getHeader();
    }
    int sampleNumber = 0;
    if (header.getSampleNames().size() > 1 || sampleName != null) {
      final int sn = header.getSampleIndex(sampleName);
      if (sn == -1) {
        if (sampleName == null) {
          throw new NoTalkbackSlimException("Need to specify a sample name for a multi-sample VCF file: " + vcf.getPath());
        } else {
          throw new NoTalkbackSlimException("Sample name not found in VCF file: " + sampleName + " : " + vcf.getPath());
        }
      }
      sampleNumber = sn;
    }
    return sampleNumber;
  }

}
