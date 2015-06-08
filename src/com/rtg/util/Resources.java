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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

/**
 * Functions for handling resources
 *
 */
public final class Resources {

  /**
   * Private to prevent instantiation.
   */
  private Resources() {
  }

  /**
   * Returns an input stream for reading the specified resource. Uses
   * the <code>ClassLoader</code> of the passed in class if available,
   * otherwise the system <code>ClassLoader</code>.
   *
   * @param c the class to get the class loader from
   * @param name the resource name
   * @return an input stream for reading the resource, or <code>null</code>
   * if the resource could not be found
   */
  public static InputStream getResourceAsStream(final Class<?> c, final String name) {
    ClassLoader cl = c.getClassLoader();
    if (cl == null) {
      cl = ClassLoader.getSystemClassLoader();
    }
    return cl.getResourceAsStream(name);
  }

  /**
   * Returns an input stream for reading the specified resource. This
   * method attempts to use the ClassLoader belonging to the Resources
   * class to read the resource. If this does not exist, then it uses
   * the system ClassLoader.
   *
   * @param name the resource name
   * @return an input stream for reading the resource, or <code>null</code>
   * if the resource could not be found
   */
  public static InputStream getResourceAsStream(final String name) {
    ClassLoader cl = Resources.class.getClassLoader();
    if (cl == null) {
      cl = ClassLoader.getSystemClassLoader();
    }
    return cl.getResourceAsStream(name);
  }

  /**
   * Returns a URL for the specified resource. This
   * method attempts to use the ClassLoader belonging to the Resources
   * class to read the resource. If this does not exist, then it uses
   * the system ClassLoader.
   *
   * @param name the resource name
   * @return a URL for the resource, or <code>null</code>
   * if the resource could not be found
   */
  public static URL getResource(final String name) {
    ClassLoader cl = Resources.class.getClassLoader();
    if (cl == null) {
      cl = ClassLoader.getSystemClassLoader();
    }
    return cl.getResource(name);
  }

  /**
   * List resources within a &quot;directory&quot;.
   * <p>
   * NOTE: output can be problematic when path exists in multiple locations.
   * @param quoteDirectoryUnquote the directory, requires trailing slash
   * @return the list of files in said directory
   * @throws URISyntaxException if it happens
   * @throws IOException if it happens
   */
  public static String[] listResources(String quoteDirectoryUnquote) throws URISyntaxException, IOException {
    final String dirWithSlash = trailingSlash(quoteDirectoryUnquote, true);
    final URL url;
    //try with and without trailing slash
    URL temp = getResource(dirWithSlash);
    if (temp == null) {
      temp = getResource(trailingSlash(quoteDirectoryUnquote, false));
    }
    url = temp;
    if (url == null) {
      return null;
    }
    if (url.getProtocol().equals("file")) {
      final ArrayList<String> strings = new ArrayList<>();
      final File dir = new File(url.toURI());
      for (final File f : dir.listFiles()) {
        strings.add(StringUtils.FS.equals("/") ? f.getPath() : f.getPath().replaceAll("\\" + StringUtils.FS, "/"));
      }
      return strings.toArray(new String[strings.size()]);
    } else if (url.getProtocol().equals("jar")) {
      final ArrayList<String> strings = new ArrayList<>();
      assert url.getPath().startsWith("file:");
      final String jarPath = url.getPath().substring(5, url.getPath().lastIndexOf('!'));
      try (JarFile jf = new JarFile(jarPath)) {
        final Enumeration<JarEntry> en = jf.entries();
        while (en.hasMoreElements()) {
          final JarEntry je = en.nextElement();
          if (je.getName().startsWith(dirWithSlash)) {
            final String rest = je.getName().substring(dirWithSlash.length());
            final int index = rest.indexOf('/');
            if (rest.length() > 0 && (index < 0 || (index > 0 && index == rest.length() - 1))) {
              strings.add(je.getName());
            }
          }
        }
      }
      return strings.toArray(new String[strings.size()]);
    }
    return null;
  }

  /**
   * @param dir directory
   * @param on true to return with a trailing slash, false for without
   * @return the directory with or without a trailing slash
   */
  static String trailingSlash(String dir, boolean on) {
    if (on) {
      if (dir.endsWith("/")) {
        return dir;
      } else {
        return dir + "/";
      }
    }
    //on == false;
    if (dir.endsWith("/")) {
      return dir.substring(0, dir.length() - 1);
    }
    return dir;
  }
}

