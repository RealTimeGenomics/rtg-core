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

package com.rtg.variant.avr;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Optional;

import com.rtg.ml.Attribute;
import com.rtg.ml.BuildClassifier;
import com.rtg.ml.BuilderFactory;
import com.rtg.ml.Dataset;
import com.rtg.ml.Instance;
import com.rtg.ml.PredictClassifier;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.StringUtils;
import com.rtg.util.ThreadAware;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.avr.AttributeExtractor.IncompatibleHeaderException;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.eval.WithInfoEvalSynchronizer;
import com.rtg.vcf.header.FormatField;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.VcfHeader;


/**
 * A machine learning model builder, bridges between the VCF world and ML package.
 */
public class MlAvrModelBuilder extends AbstractModelBuilder<MlAvrPredictModel> implements ThreadAware {

  /** Name of the configuration parameter that determines the type of ML model to build */
  public static final String PARAMETER_ML_MODEL_TYPE = "avr.subclassifier";

  /** If set, inject noise at the specified rate */
  public static final String PARAMETER_ML_INJECT_NOISE = "avr.inject-noise";

  private int mNumThreads = 1;
  private Map<String, InfoField> mCurrentInfos;
  private Map<String, FormatField> mCurrentFormats;

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

  @Override
  public void build(VcfDataset... vcfDatasets) throws IOException {

    scanVcfHeaders(vcfDatasets);
    final Annotation[] annotations = createAnnotations();
    final Attribute[] attributes = AttributeExtractor.createAttributes(annotations);
    final Dataset dataset = extractDataset(attributes, vcfDatasets);

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
    final AttributeExtractor ae = new AttributeExtractor(annotations, attributes);
    Diagnostic.info("Score for example with all values missing: " + Utils.realFormat(mlClassifier.predict(ae.getMissingValuesInstance()), 2));
    mModel = new MlAvrPredictModel(mlClassifier);
    mModel.setAttributeExtractor(ae);
    Diagnostic.userLog("Finished build");
  }

  // Create a Dataset from the pos/neg examples
  protected Dataset extractDataset(Attribute[] attributes, VcfDataset... vcfDatasets) throws IOException {
    final SimpleThreadPool stp = new SimpleThreadPool(getNumberOfThreads(), "dataset-build", false);
    final Map<VcfDataset, Dataset> datasets = new HashMap<>();
    for (final VcfDataset vcfDataset : vcfDatasets) {
      stp.execute(() -> {
        Diagnostic.userLog("Loading " + vcfDataset);
        try (final VcfReader reader = VcfReader.openVcfReader(vcfDataset.getVcfFile(), vcfDataset.ranges())) {
          final AttributeExtractor ae = new AttributeExtractor(createAnnotations(), attributes);
          try {
            ae.checkHeader(reader.getHeader());
          } catch (IncompatibleHeaderException ihe) {
            throw new NoTalkbackSlimException("The input VCF header is missing required fields:" + StringUtils.LS + ihe.getMessage());
          }
          final Dataset dataset = new Dataset(attributes);
          if (vcfDataset.classifications() == VcfDataset.Classifications.ANNOTATED) {
            if (reader.getHeader().getInfoField(WithInfoEvalSynchronizer.INFO_CALL) == null) {
              throw new NoTalkbackSlimException("The input VCF header is missing required INFO field: CALL");
            }
            while (reader.hasNext()) {
              final VcfRecord rec = reader.next();
              getClassification(rec).ifPresent(isPos -> dataset.addInstance(new Instance(ae.getInstance(rec, vcfDataset.getSampleNum()), isPos, vcfDataset.getInstanceWeight())));
            }
          } else {
            reader.forEach(rec -> dataset.addInstance(new Instance(ae.getInstance(rec, vcfDataset.getSampleNum()), vcfDataset.isPositive(), vcfDataset.getInstanceWeight())));
          }
          synchronized (datasets) {
            datasets.put(vcfDataset, dataset);
          }
          Diagnostic.userLog("Finished " + vcfDataset);
        }
      });
    }
    stp.terminate();
    Diagnostic.userLog("Merging datasets");
    final Dataset dataset = new Dataset(attributes);
    for (final VcfDataset vcfDataset : vcfDatasets) {
      dataset.addDataset(datasets.get(vcfDataset));
    }

    // Most of the time we want to make the total weight of the positive
    // instances equal to the total weight of negative instances.  We
    // only don't do it, if the user explicitly set some other weight.
    boolean reweight = true;
    for (final VcfDataset vcfDataset : vcfDatasets) {
      reweight &= vcfDataset.isReweight();
    }
    if (reweight) {
      Diagnostic.userLog("Reweighting dataset");
      dataset.reweight();
    }

    final double errorRate = Double.parseDouble(mParameters.getProperty(PARAMETER_ML_INJECT_NOISE, "-1"));
    if (errorRate > 0.0) {
      Diagnostic.userLog("Injecting noise at rate " + errorRate);
      dataset.injectMissing(errorRate);
    }

    Diagnostic.info("Total number of examples: " + dataset.size());
    Diagnostic.info("Total number of positive examples: " + dataset.totalPositives());
    Diagnostic.info("Total number of negative examples: " + dataset.totalNegatives());
    Diagnostic.info("Total weight of positive examples: " + Utils.realFormat(dataset.totalPositiveWeight(), 2));
    Diagnostic.info("Total weight of negative examples: " + Utils.realFormat(dataset.totalNegativeWeight(), 2));
    Diagnostic.info(new AttributeExtractor(createAnnotations(), attributes).missingValuesReport(dataset));
    return dataset;
  }

  /**
   * Determines classification status for a record based on vcfeval style annotation.
   *
   * @param rec the input record
   * @return true if the instance is a positive example, false for negative example, or null if no recognized status annotation is present
   */
  private Optional<Boolean> getClassification(VcfRecord rec) {
    final String status = rec.getInfo(WithInfoEvalSynchronizer.INFO_CALL);
    if (WithInfoEvalSynchronizer.STATUS_TP.equals(status)) {
      return Optional.of(true);
    } else if (WithInfoEvalSynchronizer.STATUS_FP.equals(status)) {
      return Optional.of(false);
    }
    return Optional.empty();
  }

  protected void scanVcfHeaders(VcfDataset... vcfDatasets) throws IOException {
    // Extract available annotation information from headers
    mCurrentInfos = new HashMap<>();
    mCurrentFormats = new HashMap<>();
    for (VcfDataset vcfDataset : vcfDatasets) {
      try (final VcfReader reader = VcfReader.openVcfReader(vcfDataset.getVcfFile())) {
        final VcfHeader header = reader.getHeader();
        assert header != null;
        for (final InfoField field : header.getInfoLines()) {
          final String fieldId = field.getId();
          final InfoField existing = mCurrentInfos.get(fieldId);
          // Check for header field inconsistencies (just those we think we're going to use during the build)
          if (existing != null && Arrays.binarySearch(mInfoAttributes, field.getId()) >= 0 && !field.equals(existing)) {
            throw new NoTalkbackSlimException("Info field " + fieldId + " has different definitions in different input VCFs");
          }
          mCurrentInfos.put(fieldId, field);
        }
        for (final FormatField field : header.getFormatLines()) {
          final String fieldId = field.getId();
          final FormatField existing = mCurrentFormats.get(fieldId);
          // Check for header field inconsistencies (just those we think we're going to use during the build)
          if (existing != null && Arrays.binarySearch(mFormatAttributes, field.getId()) >= 0 && !field.equals(existing)) {
            throw new NoTalkbackSlimException("Format field " + fieldId + " has different definitions in different input VCFs");
          }
          mCurrentFormats.put(field.getId(), field);
        }
      }
    }
  }

  // Create a new set of annotation according to current configuration.
  protected Annotation[] createAnnotations() {
    final List<Annotation> annotations = new ArrayList<>(mInfoAttributes.length + mFormatAttributes.length + mDerivedAttributes.length);
    if (mUseQualAttribute) {
      annotations.add(new QualAnnotation());
    }
    for (String attr : mInfoAttributes) {
      if (!mCurrentInfos.containsKey(attr)) {
        throw new NoTalkbackSlimException("Unknown INFO field: " + attr);
      }
      annotations.add(new InfoAnnotation(attr, AttributeExtractor.getCompatibleType(mCurrentInfos.get(attr).getType())));
    }
    for (String attr : mFormatAttributes) {
      if (!mCurrentFormats.containsKey(attr)) {
        throw new NoTalkbackSlimException("Unknown FORMAT field: " + attr);
      }
      annotations.add(new FormatAnnotation(attr, AttributeExtractor.getCompatibleType(mCurrentFormats.get(attr).getType())));
    }
    for (String attr : mDerivedAttributes) {
      annotations.add(new DerivedAnnotation(attr.toUpperCase(Locale.getDefault())));
    }

    return AttributeExtractor.normalizeAnnotations(annotations);
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
