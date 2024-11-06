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
  private final boolean mMakeBlacklist;
  private final boolean mInstallBlacklist;
  private final int mBlacklistThreshold;
  private final File mDirectory;
  private final int mNumberThreads;
  private final double mHashMapSizeFactor;

  /**
   * @param builder the builder
   */
  HashDistParams(HashDistParamsBuilder builder) {
    super("HashToolsParams");
    mBuildParams = builder.mBuildParams;
    mThreshold = builder.mThreshold;
    mMakeBlacklist = builder.mMakeBlacklist;
    mBlacklistThreshold = builder.mMakeBlacklist ? builder.mBlacklistThreshold : 0;
    mDirectory = builder.mDirectory;
    mNumberThreads = builder.mNumberThreads;
    mInstallBlacklist = builder.mInstallBlacklist;
    mHashMapSizeFactor = builder.mHashMapSizeFactor;
  }

  /**
   * @return true if a blacklist should be generated
   */
  public boolean makeBlacklist() {
    return mMakeBlacklist;
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
  public BuildParams buildParams() {
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
