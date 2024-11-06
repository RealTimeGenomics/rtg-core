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

package com.rtg.util;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.rtg.util.diagnostic.SlimException;

/**
 * A separate class loader implementation to allow separate copies
 * of specific classes in memory.
 *
 * This is to be used to allow JIT to optimise certain pieces of often
 * repeated code when accessing from two different locations has caused
 * optimisation to not work as well.
 */
public final class SeparateClassLoader extends ClassLoader {
  private final Set<String> mClassNames;
  private final Map<String, Class<?>> mRetrievedClasses;

  /**
   * Constructor for class loader to load separate instances of classes.
   * @param separateClasses the list of classes to use the custom class loader for.
   */
  public SeparateClassLoader(final Class<?>... separateClasses) {
    super(SeparateClassLoader.class.getClassLoader());
    mClassNames = new HashSet<>(separateClasses.length);
    for (Class<?> separateClass : separateClasses) {
      mClassNames.add(separateClass.getName());
    }
    mRetrievedClasses = new HashMap<>();
  }

  @Override
  public Class<?> loadClass(String name) throws ClassNotFoundException {
    if (mClassNames.contains(name)) {
      Class<?> clazz = mRetrievedClasses.get(name);
      if (clazz != null) {
        return clazz;
      }
      try {
        final InputStream input = Resources.getResourceAsStream(name.replace('.', '/') + ".class");
        final ByteArrayOutputStream buffer = new ByteArrayOutputStream();
        int data = input.read();

        while (data != -1) {
            buffer.write(data);
            data = input.read();
        }

        input.close();

        final byte[] classData = buffer.toByteArray();

        clazz = defineClass(name, classData, 0, classData.length);
        mRetrievedClasses.put(name, clazz);
        return clazz;
      } catch (IOException e) {
        throw new SlimException("Unable to load class: " + name);
      }
    }
    return super.loadClass(name);
  }
}
