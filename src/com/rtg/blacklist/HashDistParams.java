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
import java.io.IOException;

import com.rtg.launcher.BuildParams;
import com.rtg.launcher.ModuleParams;

/**
 * Configuration object for hash tools command
 */
public class HashDistParams extends ModuleParams {

  private final BuildParams mBuildParams;

  private final int mThreshold;
  private final boolean mBlacklist;
  private final boolean mInstallBlacklist;
  private final int mBlacklistThreshold;
  private final File mDirectory;
  private final int mNumberThreads;
  private final double mHashMapSizeFactor;

  /**
   * @param builder the builder
   */
  public HashDistParams(HashDistParamsBuilder builder) {
    super("HashToolsParams");
    mBuildParams = builder.mBuildParams;
    mThreshold = builder.mThreshold;
    mBlacklist = builder.mBlacklist;
    mBlacklistThreshold = builder.mBlacklist ? builder.mBlacklistThreshold : 0;
    mDirectory = builder.mDirectory;
    mNumberThreads = builder.mNumberThreads;
    mInstallBlacklist = builder.mInstallBlacklist;
    mHashMapSizeFactor = builder.mHashMapSizeFactor;
  }

  /**
   * @return true if a blacklist should be generated
   */
  public boolean makeBlacklist() {
    return mBlacklist;
  }

  /**
   * @return hashes with counts over this will be placed in the blacklist
   */
  public int blacklistThreshold() {
    return mBlacklistThreshold;
  }

  /**
   * @return whether generated blacklist should be installed in SDF
   */
  public boolean installBlacklist() {
    return mInstallBlacklist;
  }

  /**
   * @return the multiplier for the minimum size of the hash map
   */
  public double hashMapSizeFactor() {
    return mHashMapSizeFactor;
  }

  /**
   * @return the minimum hash count maximum. (i.e. the hash counter should be able to store counts at leas this large)
   */
  public int threshold() {
    return mThreshold;
  }

  /**
   * @return Specifies the word and step size and the sequences to generate hashes from
   */
  public BuildParams build() {
    return mBuildParams;
  }

  /**
   * @return number of threads to use
   */
  public int numberThreads() {
    return mNumberThreads;
  }

  @Override
  public File directory() {
    return mDirectory;
  }

  @Override
  public File file(String name) {
    return new File(mDirectory, name);
  }

  @Override
  public void close() throws IOException {
    mBuildParams.close();
  }

  /**
   * @return the builder
   */
  public static HashDistParamsBuilder builder() {
    return new HashDistParamsBuilder();
  }
}
