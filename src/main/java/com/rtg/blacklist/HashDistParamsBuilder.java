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
