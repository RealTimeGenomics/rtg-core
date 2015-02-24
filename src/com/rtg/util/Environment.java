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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.MissingResourceException;
import java.util.Properties;
import java.util.ResourceBundle;

import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * This class reports all the Environment variables.
 *
 */
public final class Environment {

  private static Class<?> getOperatingSystemMXBeanClass() {
    Class<?> c;
    try {
      c = Class.forName("com.sun.management.OperatingSystemMXBean");
    } catch (ClassNotFoundException e) {
      c = null;
    }
    return c;
  }
  private static final Class<?> OPERATING_SYSTEM_MX_BEAN_CLASS = getOperatingSystemMXBeanClass();

  static final String RUNTIME_FREEMEM = "runtime.freememory";
  static final String RUNTIME_MAXMEM = "runtime.maxmemory";
  static final String RUNTIME_TOTALMEM = "runtime.totalMemory";

  static final String PROCESSOR_ARCH = "processor.arch";
  static final String PROCESSOR_COUNT = "processor.count";
  static final String DEFAULT_THREADS = "runtime.defaultThreads";

  static final String OS_FREEMEM = "os.freememory";
  static final String OS_TOTALMEM = "os.totalmemory";
  static final String HOST_NAME = "host.name";

  static final String VERSION_NOT_FOUND = "<not found>";
  private Environment() { }

  /**
   * Returns HashMap with following properties
   * <ul>
   *   <li> All system property </li>
   *   <li> Detailed memory information </li>
   *   <li> Detailed OS information </li>
   * </ul>
   * This method creates the map at call time, so all values are current.
   *
   * @return <code>HashMap </code> containing all required information
   */
  public static HashMap<String, String> getEnvironmentMap() {
    final HashMap<String, String> env = new HashMap<>();
    final Properties props = System.getProperties();

    for (Map.Entry<Object, Object> entry : props.entrySet()) {
      // Get property name
      final String propName = (String) entry.getKey();
      // Get property value
      final String propValue = (String) entry.getValue();
      env.put(propName, propValue);

    }
    final Runtime rt = Runtime.getRuntime();
    env.put(RUNTIME_FREEMEM, String.valueOf(rt.freeMemory()));
    env.put(RUNTIME_MAXMEM, String.valueOf(rt.maxMemory()));
    env.put(RUNTIME_TOTALMEM, String.valueOf(rt.totalMemory()));

    final OperatingSystemMXBean bean = ManagementFactory.getOperatingSystemMXBean();
    env.put(PROCESSOR_ARCH, String.valueOf(bean.getArch()));
    env.put(PROCESSOR_COUNT, String.valueOf(getAvailableProcessors()));
    try {
      final String freeMemory = String.valueOf(getFreeMemory());
      final String totalMemory = String.valueOf(getTotalMemory());

      env.put(OS_FREEMEM, freeMemory);
      env.put(OS_TOTALMEM, totalMemory);
    } catch (final IllegalStateException e) {
      // ignore
    }
    env.put("lang.version", props.get("java.version").toString());
    env.put("vendor", props.get("java.vendor").toString());
    env.put(HOST_NAME, getHostName());
    return env;
  }

  /**
   * Method returns the free memory in system.
   *
   * @return long Free available memory in the system
   * @throws IllegalStateException if the free memory could not be determined
   */
  public static long getFreeMemory() {
    final OperatingSystemMXBean bean = ManagementFactory.getOperatingSystemMXBean();
    final Class<?> c = OPERATING_SYSTEM_MX_BEAN_CLASS == null ? bean.getClass() : OPERATING_SYSTEM_MX_BEAN_CLASS;
    try {
      final Method m = c.getMethod("getFreePhysicalMemorySize", (Class<?>[]) null);
      return (Long) m.invoke(bean, (Object[]) null);
    } catch (final IllegalAccessException | InvocationTargetException e) {
      throw new RuntimeException(e);
    } catch (final NoSuchMethodException e) {
      throw new IllegalStateException();
    }
  }


  /**
   * Returns total memory available in the system
   * @return long total memory available in the system
   * @throws IllegalStateException if the total memory could not be determined
   */
  public static long getTotalMemory() {
    final OperatingSystemMXBean bean = ManagementFactory.getOperatingSystemMXBean();
    final Class<?> c = OPERATING_SYSTEM_MX_BEAN_CLASS == null ? bean.getClass() : OPERATING_SYSTEM_MX_BEAN_CLASS;
    try {
      Method m;
      try {
       m = c.getMethod("getTotalPhysicalMemorySize", (Class<?>[]) null);
      } catch (final NoSuchMethodException e) {
       m = c.getMethod("getTotalPhysicalMemory", (Class<?>[]) null);
      }
      return (Long) m.invoke(bean, (Object[]) null);
    } catch (final IllegalAccessException | InvocationTargetException e) {
      throw new RuntimeException(e);
    } catch (final NoSuchMethodException e) {
      throw new IllegalStateException();
    }
  }

  /**
   * Returns the host name
   * @return the name
   */
  public static String getHostName() {
    try {
      return InetAddress.getLocalHost().getHostName();
    } catch (final UnknownHostException e) {
      return "<Unknown Host>";
    }
  }

  /**
   * Get number of available processors.
   * @return number of available processors.
   */
  public static int getAvailableProcessors() {
    final OperatingSystemMXBean bean = ManagementFactory.getOperatingSystemMXBean();
    if (bean != null) {
      return bean.getAvailableProcessors();
    }
    throw new IllegalStateException();
  }

  /**
   * Get the default number of threads
   * @return the number of threads defined specified by the environment if set or the number of available processors
   */
  public static int defaultThreads() {
    final String s = getEnvironmentMap().get(DEFAULT_THREADS);
    if (s == null) {
      return getAvailableProcessors();
    }
    try {
      return Integer.parseInt(s);
    } catch (NumberFormatException e) {
      throw new NoTalkbackSlimException("the system configuration for default number of threads is invalid: " + s);
    }
  }

  /**
   * Check if we have a 64-bit JVM.
   * @return true iff we have a 64 bit JVM (or at least can reliably tell if we have one).
   */
  public static boolean is64BitVM() {
    final String bits = System.getProperty("sun.arch.data.model", "?");
    return bits.equals("64");
  }

  static String getVersion(final String resource) {
    try {
      try (BufferedReader r = new BufferedReader(new InputStreamReader(Resources.getResourceAsStream(resource)))) {
        return r.readLine();
      }
    } catch (final IOException e) {
      return VERSION_NOT_FOUND;
    } catch (final NullPointerException e) {
      return VERSION_NOT_FOUND;
    }
  }

  /**
   * @return the build string
   */
  public static String getBuild() {
    return getVersion("com/rtg/slim.version");
  }

  /**
   * @return the release number as a string
   */
  public static String getRelease() {
    String version;
    final ResourceBundle environment = ResourceBundle.getBundle("com.rtg.util.Environment", Locale.getDefault());
    try {
      version = environment.getString("VERSION");
    } catch (final MissingResourceException ex) {
      version = VERSION_NOT_FOUND;
    }
    return version;
  }

  /**
   * @return the full version string
   */
  public static String getVersion() {
    return getProductName() + " / Core " + getCoreVersion();
  }

  /**
   * @return the core version string
   */
  public static String getCoreVersion() {
    final String buildtime = getVersion("com/rtg/build.time");
    final String release = getRelease();
    return release + " build " + getBuild() + " (" + buildtime + ")";
  }

  /**
   * Example: <code>"RTG Metagenomics v1.0"</code>
   * @return the product name
   */
  public static String getProductName() {
    return getVersion("com/rtg/product.name");
  }

  /**
   * Test from the command line
   * @param args command line arguments
   */
  public static void main(String... args) {
    final Map<String, String> props = Environment.getEnvironmentMap();
    for (Map.Entry<String, String> entry : props.entrySet()) {
      System.out.println(entry.getKey() + "=" + entry.getValue());
    }
  }
}

