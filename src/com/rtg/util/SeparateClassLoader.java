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
