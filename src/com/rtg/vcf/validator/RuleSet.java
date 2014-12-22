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

package com.rtg.vcf.validator;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.rtg.util.StringUtils;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfNumber;
import com.rtg.vcf.header.VcfNumberType;

/**
 * A set of simple rules that a value in the INFO or FORMAT fields must follow.
 */
public class RuleSet<T> {

  /**
   * Field types to which these rule sets can be applied.
   */
  public enum FieldType {
    /** INFO VCF field type */
    INFO,
    /** FORMAT VCF field type */
    FORMAT
  }

  /**
   * Interface to convert a String to the expected type.
   */
  public interface Converter<T> {
    /**
     * Get the value in the appropriate type.
     * @param value the value to convert from a string.
     * @return the converted value.
     */
    T getValue(String value);
  }

  /**
   * Interface to validate a value by the rule.
   */
  public interface RuleValidator<T> {
    /**
     * Check if the given value is valid for this rule.
     * @param value the value to check.
     * @return true if the value is valid, false otherwise.
     */
    boolean isValid(T value);
  }

  /**
   * Basic String implementation
   */
  public static class StringConverter implements Converter<String> {
    @Override
    public String getValue(String value) {
      return value;
    }
  }

  /**
   * Checks that a value falls within a set of known values.
   */
  private static class EnumerationRule<T> implements RuleValidator<T> {
    private final Set<T> mValidValues = new HashSet<>();
    EnumerationRule(Collection<T> values) {
      mValidValues.addAll(values);
    }
    @Override
    public boolean isValid(T value) {
      return mValidValues.contains(value);
    }
  }

  private final String mName;
  private final FieldType mType;
  private final VcfNumber mNumber;
  private final MetaType mVcfType;
  private final Converter<T> mConverter;
  protected final List<RuleValidator<T>> mRules = new ArrayList<>();

  /**
   * Constructor for a rule set.
   * @param name the name of the field.
   * @param type the type of the field.
   * @param number the VCF number object describing the number of values.
   * @param vcfType the VCF meta data type.
   * @param converter the value converter to get the String into the type required for the rule checks.
   */
  public RuleSet(String name, FieldType type, VcfNumber number, MetaType vcfType, Converter<T> converter) {
    mName = name;
    mType = type;
    mNumber = number;
    mVcfType = vcfType;
    mConverter = converter;
  }

  /**
   * Get the field name.
   * @return the field name.
   */
  public String getName() {
    return mName;
  }

  /**
   * Get the meta data type for this field.
   * @return the VCF meta data type.
   */
  public MetaType getMetaType() {
    return mVcfType;
  }

  /**
   * Get the number type for this field.
   * @return the number type.
   */
  public VcfNumber getVcfNumber() {
    return mNumber;
  }

  /**
   * Get the VCF field type for this rule.
   * @return the VCF field type.
   */
  public FieldType getVcfFieldType() {
    return mType;
  }

  /**
   * Adds an enumeration rule for this rule set.
   * @param values the string containing the comma separated values that are valid.
   * @throws RuleValidationException if there is a problem parsing the value to use.
   */
  public void addEnumerationRule(String values) throws RuleValidationException {
    final String[] vals = StringUtils.split(values, ',');
    final List<T> vs = new ArrayList<>();
    for (String val : vals) {
      vs.add(getValue(val));
    }
    mRules.add(new EnumerationRule<>(vs));
  }

  /**
   * Convert to the expected type.
   * @param value the string value to convert into the expected type.
   * @return the value in the expected type.
   * @throws RuleValidationException if the value is in the incorrect format for conversion.
   */
  protected T getValue(String value) throws RuleValidationException {
    return mConverter.getValue(value);
  }

  /**
   * Validate the record for the rules of this particular field.
   * This should only be called when the expected field is present in the record.
   * @param record the record to validate.
   * @throws RuleValidationException if record has a rule violation for the field.
   */
  public void validateRecord(VcfRecord record) throws RuleValidationException {
    if (mType == FieldType.FORMAT) {
      final List<String> sampleValues = record.getFormatAndSample().get(mName);
      for (final String sampleValue : sampleValues) {
        validateValues(StringUtils.split(sampleValue, ','), record);
      }
    } else {
      final List<String> values = record.getInfo().get(mName);
      validateValues(values.toArray(new String[values.size()]), record);
    }
  }

  private void validateValues(String[] values, VcfRecord record) throws RuleValidationException {
    if (mNumber.getNumber() == 0) {
      if (values.length != 0) {
        throw new RuleValidationException("Expected no values for the " + mType + " field " + mName + ".");
      }
    } else {
      if (values.length == 0) {
        throw new RuleValidationException("Expected values for the " + mType + " field " + mName + ".");
      }
      final int maxValues;
      final int minValues;
      if (mNumber.getNumber() > 0) {
        maxValues = mNumber.getNumber();
        minValues = maxValues;
      } else if (mNumber.getNumberType() == VcfNumberType.ALTS) {
        maxValues = record.getAltCalls().size();
        minValues = maxValues;
      } else if (mNumber.getNumberType() == VcfNumberType.REF_ALTS) {
        maxValues = record.getAltCalls().size() + 1;
        minValues = maxValues;
      } else {
        maxValues = Integer.MAX_VALUE;
        minValues = 1;
      }
      //TODO: special case for GENOTYPES number type
      if (!(values.length == 1 && VcfRecord.MISSING.equals(values[0]))) {
        if (values.length > maxValues) {
          throw new RuleValidationException("Too many values for the " + mType + " field " + mName + ".");
        }
        if (values.length < minValues) {
          throw new RuleValidationException("Too few values for the " + mType + " field " + mName + ".");
        }
      }
    }
    for (String valueString : values) {
      if (!VcfRecord.MISSING.equals(valueString)) {
        final T value = getValue(valueString);
        for (RuleValidator<T> rule : mRules) {
          if (!rule.isValid(value)) {
            //TODO: have some way of converting the rule set into a simple string for this error?
            throw new RuleValidationException("One or more values for the " + mType + " field " + mName + " is outside the expected range of values.");
          }
        }
      }
    }
  }
}
