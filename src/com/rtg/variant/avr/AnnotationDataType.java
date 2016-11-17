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

package com.rtg.variant.avr;

import com.rtg.vcf.header.MetaType;

/**
 * Data type to allow conversion between AVR types and other types.
 */
public enum AnnotationDataType {
  /**
   * Boolean/flag type.
   */
  BOOLEAN {
    @Override
    public Class<?> getClassType() {
      return Boolean.class;
    }

    @Override
    public Object stringToObjectOfType(String str) {
      return str == null ? null : Boolean.valueOf(str);
    }

    @Override
    public boolean isMetaTypeCompatible(MetaType mt) {
      return mt == MetaType.FLAG;
    }
  },
  /**
   * Integer type.
   */
  INTEGER {
    @Override
    public Class<?> getClassType() {
      return Integer.class;
    }

    @Override
    public Object stringToObjectOfType(String str) {
      return str == null ? null : Integer.valueOf(str);
    }

    @Override
    public boolean isMetaTypeCompatible(MetaType mt) {
      return mt == MetaType.INTEGER;
    }
  },
  /**
   * Double/float type.
   */
  DOUBLE {
    @Override
    public Class<?> getClassType() {
      return Double.class;
    }

    @Override
    public Object stringToObjectOfType(String str) {
      return str == null ? null : Double.valueOf(str);
    }

    @Override
    public boolean isMetaTypeCompatible(MetaType mt) {
      return mt == MetaType.FLOAT;
    }
  },
  /**
   * Generic string type.
   */
  STRING {
    @Override
    public Class<?> getClassType() {
      return String.class;
    }

    @Override
    public Object stringToObjectOfType(String str) {
      return str;
    }

    @Override
    public boolean isMetaTypeCompatible(MetaType mt) {
      return mt == MetaType.STRING || mt == MetaType.CHARACTER;
    }
  };

  /**
   * Return the java class type associated with the {@link AnnotationDataType}.
   * @return java class
   */
  public abstract Class<?> getClassType();

  /**
   * Return a type specific conversion of the given string as an object.
   * @param str string to convert
   * @return object version of string
   * @throws IllegalArgumentException if string cannot be converted
   */
  public abstract Object stringToObjectOfType(String str);

  /**
   * Test whether this {@link AnnotationDataType} is compatible with the given {@link MetaType}.
   * @param mt {@link MetaType} to test
   * @return true if types are compatible
   */
  public abstract boolean isMetaTypeCompatible(MetaType mt);
}
