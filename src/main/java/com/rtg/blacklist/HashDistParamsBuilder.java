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

package com.rtg.blacklist;

import java.io.File;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.launcher.BuildParams;

/**
 * Builder for {@link HashDistParams}
 */
@TestClass("com.rtg.blacklist.HashDistParamsTest")
public class HashDistParamsBuilder {
  BuildParams mBuildParams;
  int mThreshold;
  boolean mMakeBlacklist;
  Integer mBlacklistThreshold = 0;
  File mDirectory;
  int mNumberThreads;
  boolean mInstallBlacklist;
  double mHashMapSizeFactor = 1.0;


  /**
   * @param buildParams Specifies the word and step size and the sequences to generate hashes from
   * @return this builder for chaining purposes
   */
  public HashDistParamsBuilder buildParams(BuildParams buildParams) {
    mBuildParams = buildParams;
    return this;
  }


  /**
   * @param threshold the minimum hash count maximum. (i.e. the hash counter should be able to store counts at leas this large)
   * @return this builder for chaining purposes
   */
  public HashDistParamsBuilder threshold(int threshold) {
    mThreshold = threshold;
    return this;
  }

  /**
   * @param blacklist true if a blacklist should be generated
   * @return this builder for chaining purposes
   */
  public HashDistParamsBuilder makeBlacklist(boolean blacklist) {
    mMakeBlacklist = blacklist;
    return this;
  }

  /**
   * @param blacklistThreshold hashes with counts over this will be placed in the blacklist
   * @return this builder for chaining purposes
   */
  public HashDistParamsBuilder blacklistThreshold(Integer blacklistThreshold) {
    mBlacklistThreshold = blacklistThreshold;
    return this;
  }

  /**
   * @param directory the output directory
   * @return this builder for chaining purposes
   */
  public HashDistParamsBuilder directory(File directory) {
    mDirectory = directory;
    return this;
  }

  /**
   * @param t number of threads that should be used
   * @return this builder for chaining purposes
   */
  public HashDistParamsBuilder numberThreads(int t) {
    mNumberThreads = t;
    return this;
  }

  /**
   * @param val whether generated blacklist should be installed in SDF
   * @return this builder for chaining purposes
   */
  public HashDistParamsBuilder installBlacklist(boolean val) {
    mInstallBlacklist = val;
    return this;
  }

  /**
   * @param val the multiplier for the minimum size of the hash map
   * @return this builder for chaining purposes
   */
  public HashDistParamsBuilder hashMapSizeFactor(double val) {
    mHashMapSizeFactor = val;
    return this;
  }


  /**
   * construct the final params object from this builder
   * @return the params object
   */
  public HashDistParams create() {
    return new HashDistParams(this);
  }


  protected HashDistParamsBuilder self() {
    return this;
  }
}