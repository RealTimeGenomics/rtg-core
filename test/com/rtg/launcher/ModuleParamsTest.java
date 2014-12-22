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
package com.rtg.launcher;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.StringWriter;

import com.rtg.launcher.ModuleParams.ModuleParamsBuilder;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;


/**
 */
public class ModuleParamsTest extends TestCase {

  static final class MockModuleParamsBuilder extends ModuleParamsBuilder<MockModuleParamsBuilder> {
    @Override
    protected MockModuleParamsBuilder self() {
      return this;
    }

  }

  private static final class MockModuleParams extends ModuleParams {

    public MockModuleParams(MockModuleParamsBuilder builder) {
      super(builder);
    }
    public MockModuleParams(final CFlags flags) {
      super(flags.getName());
    }
    @Override
    public File directory() {
      return null;
    }
    @Override
    public File file(final String name) {
      return null;
    }
  }

  ModuleParams getParams(final String[] args, final String name) {
    final Appendable out = new StringWriter();
    final CFlags flags = new CFlags(name, out, null);
    flags.setFlags(args);
    final ModuleParams params = new MockModuleParams(flags);
    assertTrue(params.closed());
    return params;
  }

  public void testEquals() throws Exception {
    final ModuleParams a1 = getParams(new String[] {}, "testCliParams");
    final ModuleParams a2 = getParams(new String[] {}, "testCliParams");
    final ModuleParams b = getParams(new String[] {}, "boo");
    TestUtils.equalsHashTest(new ModuleParams[][] {{a1, a2}, {b}});

    assertEquals("testCliParams" + LS, a1.toString());
    assertEquals("testCliParams", a1.name());

    assertEquals("boo" + LS, b.toString());
    assertEquals("boo", b.name());


    a1.close();
    a2.close();
    b.close();
  }

  public void testBuilderDefaults() {
    assertEquals("ModuleParams", new MockModuleParams(new MockModuleParamsBuilder()).name());
  }

  public void testBuilder() {
    final MockModuleParamsBuilder builder = new MockModuleParamsBuilder();
    assertEquals(builder, builder.name("blah"));
  }

  public void testOmnes() {
    new TestParams(ModuleParams.class, ModuleParamsBuilder.class).check();
  }

}
