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
package com.rtg.jmx;

import static com.rtg.jmx.MonUtils.NF0;
import static com.rtg.jmx.MonUtils.NF1;
import static com.rtg.jmx.MonUtils.NF2;

import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.OperatingSystemMXBean;
import java.lang.management.RuntimeMXBean;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.rtg.util.Constants;

/**
 * Output stats from management beans.
 */
@JumbleIgnore
public class MBeanStats implements MonStats {

  private static final Class<?> OPERATING_SYSTEM_MX_BEAN_CLASS;
  static {
    Class<?> c;
    try {
      c = Class.forName("com.sun.management.OperatingSystemMXBean");
    } catch (final ClassNotFoundException e) {
      c = null;
    }
    OPERATING_SYSTEM_MX_BEAN_CLASS = c;
  }


  private static final double GB = Constants.GB;
  private static final double BILLION = 1000 * 1000 * 1000;
  private static final String LS = System.lineSeparator();

  private final MemoryMXBean mMBean;
  private final RuntimeMXBean mRBean;
  private final ThreadMXBean mTBean;
  private final OperatingSystemMXBean mOBean;

  MBeanStats() {
    mMBean = ManagementFactory.getMemoryMXBean();
    mRBean = ManagementFactory.getRuntimeMXBean();
    mTBean = ManagementFactory.getThreadMXBean();
    mOBean = ManagementFactory.getOperatingSystemMXBean();
  }

  private long getSunOsValue(String method) {
    try {
      final Method m1 = OPERATING_SYSTEM_MX_BEAN_CLASS.getMethod(method, (Class<?>[]) null);
      return (Long) m1.invoke(mOBean, (Object[]) null);
    } catch (final IllegalAccessException | InvocationTargetException e) {
      throw new RuntimeException(e);
    } catch (final NoSuchMethodException e) {
      throw new IllegalStateException();
    }
  }

  @Override
  public void addHeader(Appendable out) throws IOException {
    out.append("# Start-time   = ").append(String.valueOf(new java.util.Date(mRBean.getStartTime()))).append(LS);
    out.append("# Total-procs  = ").append(String.valueOf(mOBean.getAvailableProcessors())).append(LS);
    if (OPERATING_SYSTEM_MX_BEAN_CLASS != null) {
      out.append("# Total-mem    = ").append(NF2.format(getSunOsValue("getTotalPhysicalMemorySize") / GB)).append(" GB").append(LS);
      out.append("# Total-swap   = ").append(NF2.format(getSunOsValue("getTotalSwapSpaceSize") / GB)).append(" GB").append(LS);
    }

    out.append("# Heap-max     = ").append(NF2.format(mMBean.getHeapMemoryUsage().getMax() / GB)).append(" GB").append(LS);
    out.append("# Heap-init    = ").append(NF2.format(mMBean.getHeapMemoryUsage().getInit() / GB)).append(" GB").append(LS);
    out.append("# Nonheap-max  = ").append(NF2.format(mMBean.getNonHeapMemoryUsage().getMax() / GB)).append(" GB").append(LS);
    out.append("# Nonheap-init = ").append(NF2.format(mMBean.getNonHeapMemoryUsage().getInit() / GB)).append(" GB").append(LS);
  }

  @Override
  public void addColumnLabelsTop(Appendable out) throws IOException {
    out.append("---Up");
    if (OPERATING_SYSTEM_MX_BEAN_CLASS != null) {
      out.append(" ------OS-mem-----");
    }
    out.append(" ---Heap---- -Non-heap-- -Thrd ---OS");
    if (OPERATING_SYSTEM_MX_BEAN_CLASS != null) {
      out.append(" ---CPU");
    }
  }

  @Override
  public void addColumnLabelsBottom(Appendable out) throws IOException {
    out.append(" secs");
    if (OPERATING_SYSTEM_MX_BEAN_CLASS != null) {
      out.append("  comm   mem  swap");
    }
    out.append("  comm  used  comm  used count  load");
    if (OPERATING_SYSTEM_MX_BEAN_CLASS != null) {
      out.append("   time");
    }
  }

  @Override
  public void addColumnData(Appendable out) throws IOException {
    final int width = 5;
    MonUtils.pad(out, "" + mRBean.getUptime() / 1000, 5); // Seconds

    if (OPERATING_SYSTEM_MX_BEAN_CLASS != null) {
      out.append(" ");
      MonUtils.pad(out, NF2.format(getSunOsValue("getCommittedVirtualMemorySize") / GB), width);
      out.append(" ");
      MonUtils.pad(out, NF2.format(getSunOsValue("getFreePhysicalMemorySize") / GB), width);
      out.append(" ");
      MonUtils.pad(out, NF2.format(getSunOsValue("getFreeSwapSpaceSize") / GB), width);
    }

    out.append(" ");
    MonUtils.pad(out, NF2.format(mMBean.getHeapMemoryUsage().getCommitted() / GB), width);
    out.append(" ");
    MonUtils.pad(out, NF2.format(mMBean.getHeapMemoryUsage().getUsed() / GB), width);
    out.append(" ");
    MonUtils.pad(out, NF2.format(mMBean.getNonHeapMemoryUsage().getCommitted() / GB), width);
    out.append(" ");
    MonUtils.pad(out, NF2.format(mMBean.getNonHeapMemoryUsage().getUsed() / GB), width);
    out.append(" ");
    MonUtils.pad(out, "" + mTBean.getThreadCount(), width);
    out.append(" ");
    MonUtils.pad(out, NF1.format(mOBean.getSystemLoadAverage()), width);

    if (OPERATING_SYSTEM_MX_BEAN_CLASS != null) {
      out.append(" ");
      MonUtils.pad(out, NF0.format(getSunOsValue("getProcessCpuTime") / BILLION), width + 1);
    }
  }
}
