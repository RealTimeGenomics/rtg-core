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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import com.rtg.ml.BuildClassifier;
import com.rtg.ml.BuilderFactory;
import com.rtg.ml.Dataset;
import com.rtg.ml.Instance;
import com.rtg.ml.PredictClassifier;
import com.rtg.util.StringUtils;
import com.rtg.util.ThreadAware;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.avr.AttributeExtractor.IncompatibleHeaderException;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.header.FormatField;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfHeaderMerge;

/**
 * A machine learning model builder.
 */
public class MlAvrModelBuilder extends AbstractModelBuilder<MlAvrPredictModel> implements ThreadAware {

  /** Name of the configuration parameter that determines the type of ML model to build */
  public static final String PARAMETER_ML_MODEL_TYPE = "avr.subclassifier";

  private int mNumThreads = 1;

  /**
   * Create a Machine Learning AVR model by extracting values from the given attribute fields.
   * @param formatAttributes FORMAT field attributes to use
   * @param infoAttributes INFO field attributes to use
   * @param derivedAttributes derived field attributes to use
   */
  public MlAvrModelBuilder(String[] formatAttributes, String[] infoAttributes, String[] derivedAttributes) {
    super(formatAttributes, infoAttributes, derivedAttributes);
    mProperties.setProperty(MODEL_PROPERTY_TYPE, ModelType.ML.toString());

  }

  private double[] getMissingValuesInstance(final int size) {
    final double[] nullInstance = new double[size];
    Arrays.fill(nullInstance, Double.NaN);
    return nullInstance;
  }

  @Override
  public void build(VcfDataset... vcfDatasets) throws IOException {
    // create a merged header from all input sets
    VcfHeader mergedHeader = null;
    for (VcfDataset vcfDataset : vcfDatasets) {
      try (final VcfReader reader = VcfReader.openVcfReader(vcfDataset.getVcfFile())) {
        mergedHeader = (mergedHeader == null)
            ? reader.getHeader()
            : VcfHeaderMerge.mergeHeaders(mergedHeader, reader.getHeader(), null);
      }
    }

    // Extract available annotation information from the merged header
    final Map<String, InfoField> infos = new HashMap<>();
    final Map<String, FormatField> formats = new HashMap<>();
    assert mergedHeader != null;
    for (InfoField field : mergedHeader.getInfoLines()) {
      infos.put(field.getId(), field);
    }
    for (FormatField field : mergedHeader.getFormatLines()) {
      formats.put(field.getId(), field);
    }

    // get annotation fields
    final List<Annotation> annotations = new ArrayList<>();
    if (mUseQualAttribute) {
      annotations.add(new QualAnnotation());
    }
    for (String attr : mInfoAttributes) {
      if (!infos.containsKey(attr)) {
        throw new NoTalkbackSlimException("Unknown INFO field: " + attr);
      }
      annotations.add(new InfoAnnotation(attr, AttributeExtractor.getCompatibleType(infos.get(attr).getType())));
    }
    for (String attr : mFormatAttributes) {
      if (!formats.containsKey(attr)) {
        throw new NoTalkbackSlimException("Unknown FORMAT field: " + attr);
      }
      annotations.add(new FormatAnnotation(attr, AttributeExtractor.getCompatibleType(formats.get(attr).getType())));
    }
    for (String attr : mDerivedAttributes) {
      annotations.add(new DerivedAnnotation(attr.toUpperCase(Locale.getDefault())));
    }

    // create a Dataset from the pos/neg examples
    final AttributeExtractor ae = new AttributeExtractor(annotations.toArray(new Annotation[annotations.size()]));

    final Dataset dataset = ae.getDataset();
    boolean reweight = true;
    for (VcfDataset vcfDataset : vcfDatasets) {
      Diagnostic.userLog("Loading " + vcfDataset);
      reweight &= vcfDataset.isReweight();
      try (final VcfReader reader = VcfReader.openVcfReader(vcfDataset.getVcfFile())) {
        try {
          ae.checkHeader(reader.getHeader());
        } catch (IncompatibleHeaderException ihe) {
          throw new NoTalkbackSlimException("The input VCF header is missing required fields:" + StringUtils.LS + ihe.getMessage());
        }
        while (reader.hasNext()) {
          dataset.addInstance(new Instance(ae.getInstance(reader.next(), vcfDataset.getSampleNum()), vcfDataset.isPositive(), vcfDataset.getInstanceWeight()));
        }
      }
    }

    // Most of the time we want to make the total weight of the positive
    // instances equal to the total weight of negative instances.  We
    // only don't do it, if the user explicitly set some other weight.
    if (reweight) {
      Diagnostic.userLog("Reweighting dataset");
      dataset.reweight();
    }

    Diagnostic.info("Total number of examples: " + dataset.size());
    Diagnostic.info("Total number of positive examples: " + dataset.totalPositives());
    Diagnostic.info("Total number of negative examples: " + dataset.totalNegatives());
    Diagnostic.info("Total weight of positive examples: " + Utils.realFormat(dataset.totalPositiveWeight(), 2));
    Diagnostic.info("Total weight of negative examples: " + Utils.realFormat(dataset.totalNegativeWeight(), 2));
    Diagnostic.info(ae.missingValuesReport());

    // train ML model
    final String builderType = mParameters.getProperty(PARAMETER_ML_MODEL_TYPE, BuilderFactory.BuilderType.BAGGED.name());
    final BuildClassifier classifierBuilder;
    try {
      classifierBuilder = BuilderFactory.create(builderType, mParameters);
    } catch (IllegalArgumentException iae) {
      throw new NoTalkbackSlimException(iae.getMessage());
    }
    if (classifierBuilder instanceof ThreadAware) {
      ((ThreadAware) classifierBuilder).setNumberOfThreads(getNumberOfThreads());
    }
    Diagnostic.userLog("Starting build");
    classifierBuilder.build(dataset);

    // create predict model
    final PredictClassifier mlClassifier = classifierBuilder.getClassifier();
    Diagnostic.info("Score for example with all values missing: " + Utils.realFormat(mlClassifier.predict(getMissingValuesInstance(annotations.size())), 2));
    mModel = new MlAvrPredictModel(mlClassifier);
    mModel.setAttributeExtractor(ae);
    Diagnostic.userLog("Finished build");
  }

  @Override
  public void setNumberOfThreads(int n) {
    mNumThreads = n;
  }

  @Override
  public int getNumberOfThreads() {
    return mNumThreads;
  }
}
