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
