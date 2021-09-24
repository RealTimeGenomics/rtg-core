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

import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfNumber;

/**
 * A set of simple rules that a value in the INFO or FORMAT fields must follow.
 */
public class NumericRuleSet<T extends Number & Comparable<T>> extends RuleSet<T> {

  /**
   * Basic Integer implementation
   */
  public static class LongConverter implements Converter<Long> {
    @Override
    public Long getValue(String value) {
      return Long.valueOf(value);
    }
  }

  /**
   * Basic Double implementation
   */
  public static class DoubleConverter implements Converter<Double> {
    @Override
    public Double getValue(String value) {
      return Double.valueOf(value);
    }
  }

  /**
   * A simple no NaN checker
   */
  private static class NaNRule<T extends Number> implements RuleValidator<T> {
    @Override
    public boolean isValid(T value) {
      return !Double.isNaN(value.doubleValue());
    }
  }

  /**
   * A simple no infinity checker
   */
  private static class InfinityRule<T extends Number> implements RuleValidator<T> {
    @Override
    public boolean isValid(T value) {
      return !Double.isInfinite(value.doubleValue());
    }
  }

  /**
   * A simple greater than a given value rule.
   */
  private static class GreaterThanRule<T extends Number & Comparable<T>> implements RuleValidator<T> {
    private final T mNumber;
    /**
     * Constructor
     * @param number the number that the value must be greater than.
     */
    GreaterThanRule(T number) {
      mNumber = number;
    }
    @Override
    public boolean isValid(T value) {
      return value.compareTo(mNumber) > 0;
    }
  }

  /**
   * A simple greater than or equal to a given value rule.
   */
  private static class GreaterThanEqualRule<T extends Number & Comparable<T>> implements RuleValidator<T> {
    private final T mNumber;
    GreaterThanEqualRule(T number) {
      mNumber = number;
    }
    @Override
    public boolean isValid(T value) {
      return value.compareTo(mNumber) >= 0;
    }
  }

  /**
   * A simple less than a given value rule.
   */
  private static class LessThanRule<T extends Number & Comparable<T>> implements RuleValidator<T> {
    private final T mNumber;
    LessThanRule(T number) {
      mNumber = number;
    }
    @Override
    public boolean isValid(T value) {
      return value.compareTo(mNumber) < 0;
    }
  }

  /**
   * A simple less than or equal to a given value rule.
   */
  private static class LessThanEqualRule<T extends Number & Comparable<T>> implements RuleValidator<T> {
    private final T mNumber;
    LessThanEqualRule(T number) {
      mNumber = number;
    }
    @Override
    public boolean isValid(T value) {
      return value.compareTo(mNumber) <= 0;
    }
  }

  /**
   * Constructor for a numeric rule set.
   * @param name the name of the field.
   * @param type the type of the field.
   * @param number the VCF number object describing the number of values.
   * @param vcfType the VCF meta data type.
   * @param converter the value converter to get the String into the type required for the rule checks.
   */
  public NumericRuleSet(String name, FieldType type, VcfNumber number, MetaType vcfType, Converter<T> converter) {
    super(name, type, number, vcfType, converter);
  }

  @Override
  protected T getValue(String value) throws RuleValidationException {
    try {
      return super.getValue(value);
    } catch (NumberFormatException e) {
      throw new RuleValidationException("Value in " + getVcfFieldType() + " field " + getName() + " not in the correct number format.");
    }
  }

  /**
   * Add a no NaN rule for this rule set.
   */
  public void addNaNRule() {
    mRules.add(new NaNRule<>());
  }

  /**
   * Add a no infinity rule for this rule set.
   */
  public void addInfinityRule() {
    mRules.add(new InfinityRule<>());
  }

  /**
   * Add a greater than or equal to rule for this rule set.
   * @param value the string containing the value to check against.
   * @throws RuleValidationException if there is a problem parsing the value to use.
   */
  public void addGreaterThanEqualRule(String value) throws RuleValidationException {
    mRules.add(new GreaterThanEqualRule<>(getValue(value)));
  }

  /**
   * Add a greater than rule for this rule set.
   * @param value the string containing the value to check against.
   * @throws RuleValidationException if there is a problem parsing the value to use.
   */
  public void addGreaterThanRule(String value) throws RuleValidationException {
    mRules.add(new GreaterThanRule<>(getValue(value)));
  }

  /**
   * Add a less than or equal to rule for this rule set.
   * @param value the string containing the value to check against.
   * @throws RuleValidationException if there is a problem parsing the value to use.
   */
  public void addLessThanEqualRule(String value) throws RuleValidationException {
    mRules.add(new LessThanEqualRule<>(getValue(value)));
  }

  /**
   * Add a less than rule for this rule set.
   * @param value the string containing the value to check against.
   * @throws RuleValidationException if there is a problem parsing the value to use.
   */
  public void addLessThanRule(String value) throws RuleValidationException {
    mRules.add(new LessThanRule<>(getValue(value)));
  }
}
