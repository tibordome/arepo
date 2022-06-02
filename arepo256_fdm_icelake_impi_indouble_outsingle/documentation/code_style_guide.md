Coding Guidelines for AREPO
===========================

Code extensions
---------------

* Non-standard code extensions should always be written such that they
  can be switched off if not needed, and have no side effects on
  existing code. Normally this means that they have to be enclosed in
  conditional compilation preprocessor statements (``#ifdef``), especially
  if variables in the global structures of the code need to be
  allocated for the extension. However, if the extension's execution
  can also be controlled at run-time by simple variables, then
  consider introducing a parameterfile variable to control the
  extension. In general, the number of Makefile symbols (Config.sh
  options) to control conditional compilation should be kept to a minimum.

* Do not place any substantial piece of code belonging to your
  extension into existing functions of the code. Write your own
  functions for the code extension, and only place a function call (if
  needed bracketed by an ``#ifdef``) into the appropriate place of the
  primary code. Also, place your extension functions into separate
  source files.


General code-style principles
-----------------------------

* Code formatting: Try to be consistent with the code formatting of
  the main code, which is more or less GNU-style, with a few small
  differences. You can run the indent command with the options::

    find . -name "*.c" -o -name "*.h" | xargs clang-format -i --verbose

  on your source file(s) to make the indentation consistent.

* Name functions all in lower case as a "command" that is descriptive
  of what the function does. Words are separated by underscores, e.g.
  ``calculate_normal_vector_for_triangle(...)``.
  For all functions, input arguments should come first, output arguments
  last.

* Declare function arguments which are pointers as "const" if they are
  input-only.

* Use data types with portable size whenever needed. E.g. if you need a 32bit
  integer, use ``int32_t``; for 64bit, use ``int64_t``. (These types are part of
  standard C99 and C++11.)

* Global variables (whose use should be kept to a minimum) start with
  an upper case character. They are nouns, with words separated by
  mixed lower/upper case characters, and contain no underscores. Use
  abbreviations that are clear, e.g. ``NumForceCalculations``.

* Local variables start with lowercase, and should have descriptive
  names too, except for simple loop iterators and the like. Try to
  narrow the scope of a variable as much as possible, for example by
  declaring them inside a block where they are needed. Declaration and
  initialization should be combined in one command where possible, for
  example::

    int n = get_particle_count();

  instead of::

    int n;
    n = get_particle_count();

  The same holds for loop iteration variables. Use::

    for (int i = 0; i < 1000; ++i)

  instead of::

    int i;
    for (i = 0; i < 1000; ++i)

  If this means that you have to define the same variable several times in a
  function, that's fine.

* Code and variables that are relevant only within a single source file should
  be declared ``static`` (in C) or included in the anomymous namespace (in C++).
  This improves speed and size of the resulting code and also avoids unintended
  collisions between variable names in different compilation units.
  Most importantly, it helps the reader to distinguish between internal and
  externally visible functionality.
  (How to enforce: use function prototypes only in header files and compile with
  ``-Wstrict-prototypes``.)

* Avoid repetition of code, write a function instead. Break up long
  functions into smaller more managable pieces.

* Preprocessor macros that have arguments should be avoided whenever
  possible. If needed, their names should be fully capitalized.

* Magic numbers (including numerical constants) in the code should be
  avoided and instead be replaced by a symbolic constant declared with
  ``enum`` or ``#ifdef`` in a header file and made all uppercase.

* Address all warnings emitted by the compiler when compiling with
  ``-Wall``. Unless there are good reasons, the code should compile
  without any warning.


Code comments
-------------

Include consistent commenting in your code. The meaning of all
global variables should be commented where they are introduced,
ideally with doxygen syntax, for example::

    int MyGlobalCount;   /*!< counts the number of timesteps */

Functions should be preceeded by a brief explanation of what the
function does, including warnings/instructions, if any, about how the
function may be used. For example::

    /*!
     * Insert the point P[i] into the current tessellation. Start
     * the search at Delaunay triangle DT[t], and return the index
     * of the last accessed triangle.
     */
    int insert_point(int i, int t)
    {
      ...
    }

You do not need to state here *how* the function achieves what it
does, but this can be stated if appropriated in comments in the
function body. There, avoid superfluous comments that just reflect
what's obvious from the code anyway, like::

    do_domain_decomposition();   /* call domain decomposition */

Instead, focus on comments that help one to quickly understand/check
what the code tries to do at an algoritmic level. If complicated
formulae are implemented, try to include in a comment a reference to
the equation that is implemented.


Setting up the Eclipse IDE to work with AREPO
---------------------------------------------

Eclipse is a very powerful cross-platform development environment available for Linux, Mac, and Windows. The following instructions have been tried on a Mac, but should work also in nearly identical form on Linux.

(1) Preparation
^^^^^^^^^^^^^^^

You need to have a modern Java runtime environment on your computer, because Eclipse is written in Java::

    https://openjdk.java.net/install/

You should install the "Java SE Development Kit", which also includes the Java runtime environment.


(2) Now install Eclipse
^^^^^^^^^^^^^^^^^^^^^^^

This offers a huge number of plugins, for a large number of different languages (also including special editor extensions, for example for Markdown). These plugins can be installed from within Eclipse. What we need is "CDT", the C/C++ development platform, which can either be installed by upgrading an existing base Eclipse installation via the plug-in mechanism, or in one package together with Eclipse. I would recommend the latter approach. To do this, go to::

    https://www.eclipse.org/downloads/packages/

Select the link "Eclipse IDE for Scientific Computing" and download the build for your operating system. Then install the downloaded package.


(3) Setting up Eclipse for conveniently working with AREPO
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When starting Eclipse, you get a splash-screen offering you to look at tutorials, or carry out some basic functions, like starting a new project in various ways, etc. You can also go to the "workbench" directly (top right), which is the normal interface.

It is possible to obtain a working copy of AREPO from within Eclipse (various git functionality is built in). My preferred way is however to let Eclipse work on top of a local copy of AREPO that you obtained by cloning the repository in the command line. If you have a directory which contains a cloned copy of AREPO (which can also reside on a compute cluster if you mount the remote filesystem, if needed by tools such as sshfs which always work), then you just need to point Eclipse to it. You do this by going to "File → New → Project". Then, select "C/C++", and "Makefile Project with Existing Code". You then get the "Import Existing Code" Dialog. Type in a project name of your choice, and then select the "Existing Code Location" – here you put the directory of the cloned AREPO repository. For selecting the toolchain, I'd recommend "MacOSX GCC" on a Mac. Then click "Finish".

You should now see your new project in the left column, which you can unfold into its directory tree. Double clicking on a filename from the ``src`` directory will open the file in the editor window. You can open many files in the editor, which then appear as tabs. You can also easily have multiple editor windows, for example to see different files side by side. To this end, just try to pull windows around on the screen in the arrangement you like, they will automatically follow (or split), and you can then realize horizontal and vertical layouts of the various windows just as you find it convenient. If you want to see the same file twice (this is sometimes useful if you want to edit/see different places in the file at the same time), you can clone the corresponding editor window through "Window → Editor → Clone".

One of the main reasons for using Eclipse over a simple text-editor is that it **greatly simplifies moving aroung in large code-bases**. This can be done in many different ways. One way is to look at the "Outline" window displayed for a file (or by "unfolding" the file name in the project view), which shows the function names in the file and other prominent items. Just clicking on these jumps to the corresponding source code location. Another is to use right-click when on a function or symbol name on the editor and use functions like "Open declaration", or "References" (which shows you a list of places where this function is called, and you can directly click through those). Note that when typing new code, you can use auto-completion to automatically spell out function/variable names known in the project. Often hovering about symbols will show you helpful stuff about the symbol, for example its expected arguments, or its definition.

Besides simple searches in the local file, there are **powerful search functionalities across the whole project**. This is invoked with "Search → Search…" (or shortcut ^H, my best friend in Eclipse). In this dialog, you can do a text search (also regexp if desired) across all or selected files from the project, you can also do C/C++ syntax aware searches of various kinds. In the "File Search" dialog, you can also click a "Replace…" button at the bottom, which lets you do a search & replace in all files of a project at once. This is in principle one way to rename a function or variable through the project, but a better way is to use the tools in the "Refactor" menu, which allow you to rename functions/variables/symbols in a C/C++ syntax-aware fashion throughout the whole code-base. Take note, however, that the "Refactor" operations do not touch code which is disabled using ``#ifdef`` preprocessor directives.

To write new code and conveniently fix syntax and compile-time errors, it is **extremely useful if you can build the code from within Eclipse**. This is done with "Project → Build all" (Ctrl-B/Command-B shortcut).

To make this work, you need to have a ``Config.sh`` file in your code's directory. If none is there, you can just go to the project file list, right click on ``Template-Config.sh``, copy it, and paste it in again, changing its name to ``Config.sh`` when prompted. (If you double-click on the ".sh" file, the file opens up on a Mac in the external xcode-editor by default, but you can change this with a right click to open it in the Eclipse text editor.)

The next build error you may get is "No rule to make target `Makefile.systype'". You either need to create this file (from ``Template-Makefile.systype``) and select the system type there, or you tell Eclipse that it should use a further environment variable for the build process (``SYSTYPE="Darwin"`` for Mac). The latter can be done by going to "Project → Properties", then "C/C++ Build → Environment". There, click "Add…", and add the SYSTYPE/Darwin pair to all configurations.

The build may still fail if certain path variables (to the ``mpicc`` compiler, for example) that you set in your ``~/.bashrc`` are not yet known to Eclipse. For example, in my case, I needed to add ``/opt/local/bin/`` to ``PATH`` in Eclipse so that ``mpicc`` is found. You can do this in the same dialog as above by doing "Add…", and then using the pair, name: ``PATH``, value: ``/opt/local/bin/`` (if the option "Append variables to native environment" is set, which is default, this will be appended to your normal ``PATH`` from before ``~/.bashrc`` changes).

Now the build should work. If there are errors or warnings, you get them listed in the "Problems" window, and by clicking on them, you can easily get to the place in the source causing it. The "Console" window contains the normal output of the compiler/linker, which can also be useful to look at.

Once the build works and has finished, **the source code will be white or greyed out according to the specific configuration you have adopted in Config.sh**. If you modify this file, and rebuild the code, the selection of code adjusts accordingly. This is very useful to understand which parts of the code are actually executed for a particular configuation. A good hint is therefore to copy a ``Config.sh`` from a production setup to your source code directory, compile the code within Eclipse, and then explore it with its editor. In this way, you can understand what parts of the code are really compiled/active for this production setup.

Even if the code compiles successfully, there may be "static code analysis" warnings by Eclipse, which are labelled "Semantic Error". Many of them will be spurious in case Eclipse cannot resolve the location of certain non-standard library header files, because it doesn't know the path of these. You can recognize this by a yellow question mark on the corresponding source code line. For example, on my Mac, it wouldn't find ``mpi.h`` out of the box. So you can give Eclipse a hint about this by going to 'Project → Properties", then "C/C++ General → Code Analysis → Paths and Symbols". Then "Includes", and add a new include directory, say ``/opt/local/include/openmpi-mp/`` for the MPI header. By selecting "adding to all languages" and "all configurations" you have this generally available. After doing this, Eclipse will now know about the MPI functions in its static syntax analysis and only report real errors there. I also had to add ``/opt/local/include/`` and ``/opt/local/include/hdf5`` in this way for my installation to also properly setup HDF5, GSL, GMP, libraries, etc.

In the project explorer, when you right-click on a file, you get also access to the "Team →" functions. For example, you can simply click on a file, do "Team → Show in History", and you get the git history of the changes in this file. Clicking on one of the commits, you see what has been changed, when, and by whom, etc. Note that when you move around files in the project explorer, for example by creating new sub-directory folders in ``src/`` and moving source code files there, Eclipse will make git aware of this, i.e. this will be included in the next git commit you make for this directory.

If you mark some piece of text, you can reformat it with right-click, "Source → Format", or faster ctrl+shift+F/cmd+shift+F. But you should first set the formatting setting to our coding style! To do this, go to "Eclipse → Preferences". Then "C/C++ → Code Style → Formatter". There, select "GNU" as our base style. This at least gets close to what we want. This can be further customized in a billion ways – we should probably have someone define a default profile for this and make this available as part of AREPO's distribution.

This is all just scratching the surface about what is possible with Eclipse. It is really a huge help in working with a large code like AREPO. Once you get used to it you will wonder how you ever managed to work just with a plain text editor, so I think it is really worthwhile to get to know it.

