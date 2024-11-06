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

import java.lang.reflect.InvocationTargetException;

import junit.framework.TestCase;

/**
 * Test for loader of separate classes.
 */
public class SeparateClassLoaderTest extends TestCase {

  /**
   * Dummy Interface
   */
  public interface IFace {
    String methodA();
    String methodB();
  }

  /**
   * Dummy Abstract Class
   */
  public abstract static class AbstractClass implements IFace {
    @Override
    public String methodA() {
      return methodC();
    }
    abstract String methodC();
  }

  /**
   * Dummy Implementation Class
   */
  public static final class ImplClass extends AbstractClass {
    @Override
    public String methodB() {
      return methodC();
    }
    @Override
    String methodC() {
      return "foooooo";
    }
  }

  public void test() throws InstantiationException, IllegalAccessException, ClassNotFoundException, IllegalArgumentException, InvocationTargetException, NoSuchMethodException, SecurityException {
    final IFace normal = new ImplClass();
    assertTrue(normal instanceof AbstractClass);
    assertTrue(normal instanceof ImplClass);
    final SeparateClassLoader loader = new SeparateClassLoader(ImplClass.class, AbstractClass.class);
    final IFace separate = (IFace) loader.loadClass(ImplClass.class.getName()).getConstructor().newInstance();
    assertFalse(separate instanceof AbstractClass);
    assertFalse(separate instanceof ImplClass);
    final IFace separate2 = (IFace) loader.loadClass(ImplClass.class.getName()).getConstructor().newInstance();
    assertFalse(separate2 instanceof AbstractClass);
    assertFalse(separate2 instanceof ImplClass);
    assertEquals("foooooo", normal.methodA());
    assertEquals("foooooo", normal.methodB());
    assertEquals("foooooo", separate.methodA());
    assertEquals("foooooo", separate.methodB());
    assertEquals("foooooo", separate2.methodA());
    assertEquals("foooooo", separate2.methodB());
  }
}
