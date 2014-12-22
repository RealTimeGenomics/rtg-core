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
package com.rtg.report;

import java.io.IOException;
import java.io.StringWriter;
import java.util.HashMap;
import java.util.Map;

import org.apache.velocity.VelocityContext;
import org.apache.velocity.app.Velocity;
import org.apache.velocity.app.VelocityEngine;
import org.apache.velocity.runtime.resource.loader.ClasspathResourceLoader;

import com.rtg.util.HtmlReportHelper;

/**
 * tools for template reports using Velocity engine
 */
public final class VelocityReportUtils {
  private VelocityReportUtils() { }

  private static final String DEFAULT_TEMPLATE = "default.vm";

  static final String TEMPLATE_DIR = "com/rtg/report/resources";
  static {
    Velocity.setProperty(VelocityEngine.RESOURCE_LOADER, "class");
    Velocity.setProperty("class.resource.loader.class", ClasspathResourceLoader.class.getCanonicalName());
    Velocity.setProperty("class.resource.loader.path", "/");
    Velocity.setProperty("runtime.log.logsystem.class", RtgVelocityLogChute.class.getCanonicalName());
    Velocity.init();
  }

  /**
   * Wrap the default RTG template around some body text.
   * @param bodyText the body text
   * @param title title of page/report
   * @param hrh a HTML resources helper
   * @return the processed report
   * @throws IOException if an IO error occurs
   */
  public static String wrapDefaultTemplate(String bodyText, String title, HtmlReportHelper hrh) throws IOException {
    final Map<String, String> data = new HashMap<>();
    data.put("body", bodyText);
    data.put("title", title);
    data.put("resourceDir", hrh.getResourcesDirName());
    hrh.copyResources(TEMPLATE_DIR + "/rtg.css", TEMPLATE_DIR + "/rtg_logo.png", TEMPLATE_DIR + "/table.css");
    return processTemplate(DEFAULT_TEMPLATE, data);
  }

  /**
   * process a template, producing final output
   * @param template template for report
   * @param replacements data for report
   * @return the processed report
   * @throws IOException if an IO error occurs
   */
  public static String processTemplate(String template, Map<String, ?> replacements) throws IOException {
    final VelocityContext vc = new VelocityContext();
    for (Map.Entry<String, ?> me : replacements.entrySet()) {
      vc.put(me.getKey(), me.getValue());
    }
    final String ret;
    try (final StringWriter sw = new StringWriter()) {
      Velocity.getTemplate(TEMPLATE_DIR + "/" + template).merge(vc, sw);
      sw.flush();
      ret = sw.toString();
    }
    return ret;
  }
}
