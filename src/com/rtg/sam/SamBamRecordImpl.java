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
package com.rtg.sam;


/**
 * This gives access to all the fields in a SAM/BAM record.
 * It defines named getters for most of the fields and flags.
 *
 */
public class SamBamRecordImpl implements SamBamRecord {

  protected SamBamReader mReader;

  SamBamRecordImpl(SamBamReader source) {
    mReader = source;
  }

  /**
   * Returns the read name as a string.
   * @return the read name
   */
  @Override
  public String getReadName() {
    return mReader.getField(SamBamConstants.QNAME_FIELD);
  }

  /**
   * This is equivalent to getting the read name and parsing it as an
   * integer, but much more efficient.
   * It should be called only when you are sure the read name
   * will be a valid integer.
   *
   * @return The read name, converted to an integer.
   */
  @Override
  public int getReadNumber() {
    return mReader.getIntField(SamBamConstants.QNAME_FIELD);
  }

  /**
   * This returns all the flags as bits within an int.
   * The SAM_READ_IS_... constants can be used to test the individual
   * flags.
   *
   * @return All the flags.
   */
  @Override
  public int getFlags() {
    return mReader.getIntField(SamBamConstants.FLAG_FIELD);
  }

  /**
   * Same as in {@link SamBamReader#hasAttribute(String) SamBamReader}.
   *
   * @param tag The two-character tag of the attribute you are looking for.
   * @return true if the current record has that attribute.
   */
  @Override
  public boolean hasAttribute(String tag) {
    return mReader.hasAttribute(tag);
  }

  /**
   * @see SamBamReader#getAttributeType(String)
   *
   * @param tag The two-character tag you are looking for.
   * @return the type:value string associated with that tag, or null.
   */
  @Override
  public char getAttributeType(String tag) {
    return mReader.getAttributeType(tag);
  }

  /**
   * @see SamBamReader#getFieldNumFromTag(java.lang.String)
   *
   * @param tag the tag
   * @return the field num
   */
  @Override
  public int getFieldNumFromTag(String tag) {
    return mReader.getFieldNumFromTag(tag);
  }

  @Override
  public String[] getAttributeTags() {
    return mReader.getAttributeTags();
  }

  /**
   * This looks for an optional attribute with the given <code>tag</code>
   * and returns its value as a string/integer/double object, or null if
   * the alignment record does not contain the requested attribute.
   *
   * Pre: there must be a current alignment record.
   *
   * @param tag The two-character tag you are looking for.
   * @return the value associated with that tag, or null.
   */
  @Override
  public Object getAttributeValue(String tag) {
    return mReader.getAttributeValue(tag);
  }

  /**
   * This looks for an optional attribute with the given <code>tag</code>
   * and returns the value part as an integer.
   * It returns <code>Integer.MIN_VALUE</code>
   * if the current alignment record does not contain the requested attribute.
   *
   * Pre: there must be a current alignment record.
   *
   * @param tag The two-character tag you are looking for.
   * @return the type:value string associated with that tag, or null.
   */
  @Override
  public int getIntAttribute(String tag) {
    return mReader.getIntAttribute(tag);
  }

  /**
   * @return the number of fields in the current alignment record
   */
  @Override
  public int getNumFields() {
    return mReader.getNumFields();
  }

  /**
   * Get a string version of one field of the current alignment record.
   * Pre: <code>fieldNum</code> is at least 0 and is less than <code>getNumFields()</code>.
   *
   * @param fieldNum the zero-based number of a field in the current record
   * @return a non-null String.
   */
  @Override
  public String getField(int fieldNum) {
    return mReader.getField(fieldNum);
  }

  @Override
  public byte[] getFieldBytes(int fieldNum) {
    return mReader.getFieldBytes(fieldNum);
  }



  /**
   * Get an integer version of one field of the current alignment record.
   * Pre: <code>fieldNum</code> is at least 0 and is less than <code>getNumFields()</code>.
   *
   * @param fieldNum the zero-based number of a field in the current record
   * @return an integer value, or -1 if that field contains a non-integer.
   */
  @Override
  public int getIntField(int fieldNum) {
    return mReader.getIntField(fieldNum);
  }
}
