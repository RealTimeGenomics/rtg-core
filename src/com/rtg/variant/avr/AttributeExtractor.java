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

import static com.rtg.util.StringUtils.LS;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeSet;

import com.rtg.ml.Attribute;
import com.rtg.ml.Dataset;
import com.rtg.ml.MlDataType;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.LineWriter;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.annotation.AnnotationDataType;
import com.rtg.vcf.annotation.DerivedAnnotations;
import com.rtg.vcf.header.FormatField;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;

/**
 * Provides conversion from VCF records into instances for the machine learning.
 * Canonicalization of annotations/attributes to ensure <code>Object[]</code> has consistent order (N elements).
 */
public class AttributeExtractor {

  /**
   */
  public static final class IncompatibleHeaderException extends Exception {
    /**
     * Construct an exception with given message.
     * @param message exception message.
     */
    public IncompatibleHeaderException(String message) {
      super(message);
    }
  }

  private final Annotation[] mAnnotations;
  private final transient long[] mMissingValueCounts;
  private final Attribute[] mAttributes;

  private AttributeExtractor(int numAnnotations) {
    mAnnotations = new Annotation[numAnnotations];
    mMissingValueCounts = new long[numAnnotations];
    mAttributes = new Attribute[numAnnotations];
  }

  /**
   * Constructs an attribute extractor for the given set of annotations.
   * @param annotations annotations to process
   */
  public AttributeExtractor(Annotation... annotations) {
    final TreeSet<Annotation> set = new TreeSet<>(new Comparator<Annotation>() {
      @Override
      public int compare(Annotation o1, Annotation o2) {
        return o1.getName().compareTo(o2.getName());
      }
    });
    for (Annotation att : annotations) {
      if (att == null) {
        throw new NullPointerException("null attribute given");
      }
      if (set.contains(att)) {
        throw new IllegalArgumentException("Duplicate attribute: " + att.getName());
      }
      set.add(att);
    }

    mAnnotations = set.toArray(new Annotation[set.size()]);
    mMissingValueCounts = new long[mAnnotations.length];
    mAttributes = new Attribute[mAnnotations.length];
    initAttributes();
  }

  private void initAttributes() {
    for (int i = 0; i < mAttributes.length; i++) {
      mAttributes[i] = new Attribute(mAnnotations[i].getName(), getMlDataType(mAnnotations[i].getType()));
    }
  }

  /**
   * Check that the annotations are compatible with the VCF header.
   *
   * @param header a VCF header
   * @throws IncompatibleHeaderException if the head and attributes don't match
   */
  public void checkHeader(VcfHeader header) throws IncompatibleHeaderException {
    final StringBuilder exceptionMessages = new StringBuilder();
    for (Annotation anno : mAnnotations) {
      final String message = anno.checkHeader(header);
      if (message != null) {
        exceptionMessages.append(message).append(LS);
      }
    }
    if (exceptionMessages.length() != 0) {
      throw new IncompatibleHeaderException(exceptionMessages.toString());
    }
  }

  /*
  protected static boolean isTypeCompatible(MetaType mt, AvrDataType type) {
    switch (mt) {
      case INTEGER:
        return type == AvrDataType.INTEGER;
      case FLOAT:
        return type == Avr;
      case STRING:
        return type == String.class || (type != null && type.isEnum());
      case CHARACTER:
        return type == String.class;
      case FLAG:
        return type == Boolean.class;
      default:
        return false;
    }
  }
  */

  /**
   * Returns empty dataset that complies with the attributes.
   * @return a dataset with attributes set up.
   */
  Dataset getDataset() {

    return new Dataset(mAttributes);
  }

  protected static MlDataType getMlDataType(AnnotationDataType dataType) {
    switch (dataType) {
      case BOOLEAN:
        return MlDataType.BOOLEAN;
      case INTEGER:
        return MlDataType.INTEGER;
      case DOUBLE:
        return MlDataType.DOUBLE;
      case STRING:
        return MlDataType.STRING;
      default:
        throw new IllegalArgumentException("Unrecognised Data Type: " + dataType);
    }
  }

  /**
   * Returns an instance object array for the given VCF record.
   * @param record a VCF record
   * @param sampleNumber the sample to extract value from.
   * @return array of attribute values
   */
  public double[] getInstance(VcfRecord record, int sampleNumber) {
    final double[] res = new double[mAnnotations.length];
    for (int i = 0; i < res.length; i++) {
      try {
        res[i] = mAttributes[i].encodeValue(mAnnotations[i].getValue(record, sampleNumber));
      } catch (final NumberFormatException e) {
        throw new NoTalkbackSlimException("Problem parsing a number in a VCF record:\n" + record.toString() + "\n" + e.getMessage());
      }
      if (Double.isNaN(res[i])) {
        mMissingValueCounts[i]++;
      }
    }
    return res;
  }

  /**
   * Return a summary of the number of missing values.
   * @return missing values summary
   */
  public String missingValuesReport() {
    final String[] nums = new String[mMissingValueCounts.length];
    int maxLengthNum = 0;
    int maxLengthName = 0;
    for (int k = 0; k < nums.length; k++) {
      nums[k] = String.valueOf(mMissingValueCounts[k]);
      if (nums[k].length() > maxLengthNum) {
        maxLengthNum = nums[k].length();
      }
      final int nameLength = mAnnotations[k].getName().length();
      if (nameLength > maxLengthName) {
        maxLengthName = nameLength;
      }
    }
    final StringBuilder sb = new StringBuilder("Number of examples with missing values:").append(LS);
    for (int k = 0; k < mAnnotations.length; k++) {
      final Annotation att = mAnnotations[k];
      final String name = att.getName();
      sb.append("  ").append(name);
      for (int j = name.length(); j <= maxLengthName; j++) {
        sb.append(" ");
      }
      for (int j = nums[k].length(); j <= maxLengthNum; j++) {
        sb.append(" ");
      }
      sb.append(nums[k]).append(LS);
    }
    return sb.toString();
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    for (Annotation att : mAnnotations) {
      if (sb.length() != 0) {
        sb.append(",");
      }
      sb.append(att.getName()).append("(").append(att.getType()).append(")");
    }
    return sb.toString();
  }

  /**
   * Creates an {@link AttributeExtractor} form the given input stream.
   * @param is an input stream
   * @return an {@link AttributeExtractor} from stream contents
   * @throws IOException if an error occurs reading data
   */
  public static AttributeExtractor load(InputStream is) throws IOException {
    final DataInputStream dis = new DataInputStream(is);
    final int numAnnotations = dis.readInt();
    final AttributeExtractor ae = new AttributeExtractor(numAnnotations);
    for (int i = 0; i < numAnnotations; i++) {
      ae.mAnnotations[i] = AnnotationLoader.load(dis);
    }
    ae.initAttributes();
    return ae;
  }

  /**
   * Saves this {@link AttributeExtractor} to the given output stream.
   * @param os output stream to write to
   * @throws IOException if an error occurs writing data
   */
  public void save(OutputStream os) throws IOException {
    final DataOutputStream dos = new DataOutputStream(os);
    dos.writeInt(mAnnotations.length);
    try {
      for (Annotation anno : mAnnotations) {
        anno.save(dos);
      }
    } finally {
      dos.flush();
    }
  }

  /**
   * @param args command line arguments
   * @throws Exception when an exception occurs
   */
  public static void main(String[] args) throws Exception {
    final CFlags flags = new CFlags();
    CommonFlagCategories.setCategories(flags);
    flags.setDescription("Generate");
    flags.registerRequired('i', "input", File.class, "FILE", "input VCF file to read from").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerRequired('o', "output", File.class, "FILE", "output ARFF file").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional("info-annotations", String.class, "STRING", "comma separated list of info attributes").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional("format-annotations", String.class, "STRING", "comma separated list of format attributes").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional("sample", String.class, "STRING", "the name of the sample to select (required when using multi-sample VCF files)").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    final List<String> derivedRange = new ArrayList<>();
    for (final DerivedAnnotations derived : DerivedAnnotations.singleValueAnnotations()) {
      derivedRange.add(derived.toString());
    }
    flags.registerOptional("derived-annotations", String.class, "STRING", "derived fields to use in model").setParameterRange(derivedRange).setRangeList(true).setCategory(CommonFlagCategories.REPORTING);
    if (!flags.setFlags(args)) {
      return;
    }

    final File vcfFile = (File) flags.getValue("input");
    final File output = (File) flags.getValue("output");

    final String[] infoAnnos = flags.isSet("info-annotations") ? ((String) flags.getValue("info-annotations")).split(",") : new String[0];
    final String[] formatAnnos = flags.isSet("format-annotations") ? ((String) flags.getValue("format-annotations")).split(",") : new String[0];
    final String[] derivedAnnos;
    if (flags.isSet("derived-annotations")) {
      final List<String> derived = new ArrayList<>();
      for (final Object field : flags.getValues("derived-annotations")) {
        derived.add(field.toString());
      }
      derivedAnnos = derived.toArray(new String[derived.size()]);
    } else {
      derivedAnnos = new String[0];
    }

    final String sampleName = (String) flags.getValue("sample");

    try {
      try (final VcfReader reader = VcfReader.openVcfReader(vcfFile)) {
        final List<Annotation> annotations = new ArrayList<>();

        final VcfHeader header = reader.getHeader();

        int sampleNumber = 0;
        if (header.getSampleNames().size() > 1 || sampleName != null) {
          final Integer sn = header.getSampleIndex(sampleName);
          if (sn == null) {
            if (sampleName == null) {
              throw new NoTalkbackSlimException("Need to specify a sample name for a multi-sample VCF file.");
            } else {
              throw new NoTalkbackSlimException("Sample name not found in VCF file: " + sampleName);
            }
          }
          sampleNumber = sn;
        }

        final Map<String, InfoField> infos = new HashMap<>();
        for (InfoField field : header.getInfoLines()) {
          infos.put(field.getId(), field);
        }
        final Map<String, FormatField> formats = new HashMap<>();
        for (FormatField field : header.getFormatLines()) {
          formats.put(field.getId(), field);
        }

        annotations.add(new QualAnnotation());
        for (String anno : infoAnnos) {
          final InfoField field = infos.get(anno);
          if (field == null) {
            throw new IncompatibleHeaderException("INFO annotation not in VCF file: " + anno);
          }
          annotations.add(new InfoAnnotation(anno, getCompatibleType(field.getType())));
        }
        for (String anno : formatAnnos) {
          final FormatField field = formats.get(anno);
          if (field == null) {
            throw new IncompatibleHeaderException("FORMAT annotation not in VCF file: " + anno);
          }
          annotations.add(new FormatAnnotation(anno, getCompatibleType(field.getType())));
        }
        for (String anno : derivedAnnos) {
          annotations.add(new DerivedAnnotation(anno.toUpperCase(Locale.getDefault())));
        }

        final AttributeExtractor ae = new AttributeExtractor(annotations.toArray(new Annotation[annotations.size()]));
        ae.checkHeader(header);

        try (final LineWriter arffWriter = new LineWriter(new FileWriter(output))) {
          arffWriter.writeln("@relation '" + vcfFile.getPath() + "'");
          arffWriter.writeln("");

          for (Annotation anno : annotations) {
            arffWriter.writeln(getArffAttributeText(anno));
          }

          arffWriter.writeln("");
          arffWriter.writeln("@data");
          arffWriter.writeln("");

          while (reader.hasNext()) {
            final double[] instance = ae.getInstance(reader.next(), sampleNumber);
            for (int i = 0; i < instance.length; i++) {
              if (i != 0) {
                arffWriter.write(",");
              }
              arffWriter.write(getValueAsArffString(instance[i]));
            }
            arffWriter.writeln("");
          }
          arffWriter.flush();
        }
      }
    } catch (IncompatibleHeaderException ihe) {
      System.err.println(ihe.getMessage());
    }
  }

  protected static AnnotationDataType getCompatibleType(MetaType mt) {
    switch (mt) {
      case INTEGER:
        return AnnotationDataType.INTEGER;
      case FLOAT:
        return AnnotationDataType.DOUBLE;
      case STRING:
      case CHARACTER:
        return AnnotationDataType.STRING;
      case FLAG:
        return AnnotationDataType.BOOLEAN;
      default:
        return null;
    }
  }

  protected static String getArffAttributeText(Annotation annotation) {
    final StringBuilder sb = new StringBuilder("@attribute ");
    sb.append(annotation.getName());
    if (annotation.getType() == AnnotationDataType.INTEGER
        || annotation.getType() == AnnotationDataType.DOUBLE) {
      sb.append(" numeric");
    } else if (annotation.getType() == AnnotationDataType.BOOLEAN) {
      sb.append(" {").append(Boolean.TRUE.toString()).append(",").append(Boolean.FALSE.toString()).append("}");
    } else {
      sb.append(" string");
    }
    return sb.toString();
  }

  protected static String getValueAsArffString(Object value) {
    return value == null ? "?" : value.toString(); //.replace(',', '_');
  }
}
