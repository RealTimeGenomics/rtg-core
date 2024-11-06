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
package com.rtg.index.hash.ngs;


import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

import com.rtg.index.hash.ngs.instances.CG2Maska1;
import com.rtg.index.hash.ngs.instances.CG2Maska11;
import com.rtg.index.hash.ngs.instances.CG2Maska15;
import com.rtg.index.hash.ngs.instances.CG2Maskw18;
import com.rtg.index.hash.ngs.instances.CG2Maskw18a1;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Provides access to the set of mask factory implementations.
 *
 */
public final class FactoryUtil {

  private FactoryUtil() { } //prevent instantiation.

  private static final Map<String, HashFunctionFactory> MASK_FACTORIES = new HashMap<>();

  static {
    putFactory("CGMaska0", com.rtg.index.hash.ngs.instances.CGMaska0.FACTORY);
    putFactory("CGMaska1b1", com.rtg.index.hash.ngs.instances.CGMaska1b1.FACTORY);
    putFactory("CGMaska15b1", com.rtg.index.hash.ngs.instances.CGMaska15b1.FACTORY);
    putFactory("CGMaska1b1alt", com.rtg.index.hash.ngs.instances.CGMaska1b1alt.FACTORY);
    putFactory("CGMaska15b1alt", com.rtg.index.hash.ngs.instances.CGMaska15b1alt.FACTORY);
    putFactory("CG2Maska1", CG2Maska1.FACTORY);
    putFactory("CG2Maska11", CG2Maska11.FACTORY);
    putFactory("CG2Maska15", CG2Maska15.FACTORY);
    putFactory("CG2Maskw18", CG2Maskw18.FACTORY);
    putFactory("CG2Maskw18a1", CG2Maskw18a1.FACTORY);
    //-------------- public options
    putFactory("cg1", com.rtg.index.hash.ngs.instances.CGMaska15b1.FACTORY);
    putFactory("cg1-fast", com.rtg.index.hash.ngs.instances.CGMaska1b1.FACTORY);
    putFactory("cg2", CG2Maskw18a1.FACTORY);
  }

  private static void putFactory(String name, HashFunctionFactory factory) {
    MASK_FACTORIES.put(name.toLowerCase(Locale.getDefault()), factory);
  }

  /**
   * @param mask Class name of the mask class.
   * @return check if the mask exists
   */
  public static boolean checkMaskExists(final String mask) {
    if (mask != null) {
      return MASK_FACTORIES.containsKey(mask.toLowerCase(Locale.getDefault()));
    }
    return false;
  }

  /**
   * @param mask Class name of the mask class.
   * @return a factory which will create a new a <code>HashFunction</code>.
   */
  public static HashFunctionFactory hashFunction(final String mask) {
    final String m = mask == null ? null : mask.toLowerCase(Locale.getDefault());
    final HashFunctionFactory factory = MASK_FACTORIES.get(m);
    if (factory == null) {
      throw new NoTalkbackSlimException(ErrorType.INVALID_MASK, mask);
    }
    return factory;
  }

}

