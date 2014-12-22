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

package com.rtg.relation;

import java.util.Arrays;
import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.Utils;

/**
 * Defines a single relationship between two genomes
 */
@TestClass("com.rtg.relation.GenomeRelationshipsTest")
public class Relationship {

  /**
   * Enumeration for relationship types. Add extras as needed.
   */
  public enum RelationshipType {
    /** The first genome is a parent of the second genome. */
    PARENT_CHILD,
    /** The second genome is derived from the first genome (e.g. normal / tumor, or cell line). */
    ORIGINAL_DERIVED,
  }

  private final String mGenome1;
  private final String mGenome2;
  private final RelationshipType mType;
  private final Map<String, String> mProperties;

  /**
   * @param genome1 first genome in relationship
   * @param genome2 second genome in relationship
   * @param type type of relationship
   */
  public Relationship(String genome1, String genome2, RelationshipType type) {
    mGenome1 = genome1;
    mGenome2 = genome2;
    mType = type;
    mProperties = new TreeMap<>();
  }

  /**
   * Get genome 1.
   * @return Returns the first genome in the relationship.
   */
  public String first() {
    return mGenome1;
  }

  /**
   * Get genome 2.
   * @return Returns the second genome in the relationship.
   */
  public String second() {
    return mGenome2;
  }

  /**
   * @return Returns type of relationship
   */
  public RelationshipType type() {
    return mType;
  }

  /**
   * Sets an arbitrary property for this relationship
   * @param propname name of property
   * @param val value for property
   */
  public void setProperty(String propname, String val) {
    mProperties.put(propname, val);
  }

  /**
   * query for a property
   * @param propname name of property
   * @return value of property
   */
  public String getProperty(String propname) {
    return mProperties.get(propname);
  }

  /**
   * Get the contamination level for this relationship or null if no contamination
   * is specified.
   *
   * @return contamination level
   */
  public Double getContamination() {
    final String contamStr = getProperty("contamination");
    if (contamStr == null) {
      return null;
    }
    return Double.parseDouble(contamStr);
  }

  @Override
  public boolean equals(Object o) {
    if (!(o instanceof Relationship)) {
      return false;
    }
    final Relationship p = (Relationship) o;
    final boolean propEquals;
    propEquals = mProperties.equals(p.mProperties);
    return mGenome1.equals(p.mGenome1) && mGenome2.equals(p.mGenome2) && mType.ordinal() == p.mType.ordinal() && propEquals;
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(mProperties.hashCode(), Utils.pairHash(mType.ordinal(), Utils.pairHash(mGenome1.hashCode(), mGenome2.hashCode())));
  }

  @Override
  public String toString() {
    return mType.toString() + " (" + mGenome1 + "-" + mGenome2 + ")" + (!mProperties.isEmpty() ? (" :: " + mProperties.toString()) : "");
  }


  /**
   * Things which may accept or reject a Relationship
   */
  public interface RelationshipFilter {
    /**
     * Returns true if the relationship is accepted
     * @param relationship the relationship to test
     * @return true if the relationship is accepted
     */
    boolean accept(Relationship relationship);
  }

  /** Accepts relationships that are rejected by a sub-filter */
  public static class NotFilter implements RelationshipFilter {
    private final RelationshipFilter mFilter;

    /**
     * Accepts relationships of a certain type
     * @param filter the type to match
     */
    public NotFilter(RelationshipFilter filter) {
      mFilter = filter;
    }

    @Override
    public boolean accept(Relationship relationship) {
      return !mFilter.accept(relationship);
    }
  }

  /** Accepts relationships of a certain type */
  public static class RelationshipTypeFilter implements RelationshipFilter {
    private final RelationshipType mType;

    /**
     * Accepts relationships of a certain type
     * @param type the type to match
     */
    public RelationshipTypeFilter(RelationshipType type) {
      mType = type;
    }

    @Override
    public boolean accept(Relationship relationship) {
      return mType == relationship.type();
    }
  }

  /** Accepts relationships where a specified genome is in the first position */
  public static class FirstInRelationshipFilter implements RelationshipFilter {
    private final String mGenome;

    /**
     * Accepts relationships with a genome in the first position
     * @param genome the genome that must be in the first position
     */
    public FirstInRelationshipFilter(String genome) {
      mGenome = genome;
    }

    @Override
    public boolean accept(Relationship relationship) {
      return mGenome.equals(relationship.first());
    }
  }

  /** Accepts relationships where a specified genome is in the second position */
  public static class SecondInRelationshipFilter implements RelationshipFilter {
    private final String mGenome;

    /**
     * Accepts relationships with a genome in the second position
     * @param genome the genome that must be in the second position
     */
    public SecondInRelationshipFilter(String genome) {
      mGenome = genome;
    }

    @Override
    public boolean accept(Relationship relationship) {
      return mGenome.equals(relationship.second());
    }
  }

  /**
   * Filters relationships based on given samples. Only relationships where both genomes match are accepted.
   */
  public static class SampleRelationshipFilter implements RelationshipFilter {

    private final Collection<String> mSamples;

    /**
     * @param samples samples to keep relationships for
     */
    public SampleRelationshipFilter(String... samples) {
      mSamples = Arrays.asList(samples);
    }

    /**
     * @param samples samples to keep relationships for
     */
    public SampleRelationshipFilter(Set<String> samples) {
      mSamples = samples;
    }

    @Override
    public boolean accept(Relationship relationship) {
      return mSamples.contains(relationship.first()) && mSamples.contains(relationship.second());
    }
  }
}
