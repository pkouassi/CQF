%% usmthesis.cls 2016/04/13 version V1.66
%%
%% NOTE: This text file uses UNIX line feed conventions. When (human)
%% reading this file on other platforms, you may have to use a text
%% editor that can handle lines terminated by the UNIX line feed
%% character (0x0A).
%%
%%
%% This is the unofficial USM thesis LaTeX class for postgraduate research
%% students of Universiti Sains Malaysia, Penang, Malaysia, formatted
%% according to the Guide to the Preparation, Submission and Examination
%% of Theses, published by IPS USM.  It supercedes usmthesis.sty: much of
%% the code is essentially the same, although usmthesis is now a class
%% extending from the book class instead of a style.  usmthesis also provides
%% some class options for font choice for body text and url, as well as
%% boxed (or not) figures.
%%
%%
%% The latest version and documentation of usmthesis can be obtained at
%% http://liantze.penguinattack.org/latextypesetting.html#usmthesis.
%% contact liantze AT gmail DOT com if you have any questions or
%% assistance or feedback.  Apologies if I can't respond immediately, but
%% I'll always try.
%%
%%
%% Copyright (c) 2005-2016 by Lim Lian Tze (OH MY GOODNESS HAS IT BEEN THAT LONG)
%%
%% Current maintainer: Lim Lian Tze (liantze@gmail.com)
%
%  BRIEF CHANGELOG
%  ===============
% v1.65 Apr 13 2016
%     - First pages of appendix chapters don't have footers now.
%
% v1.65 Mar 22 2016
%     - Titles on cover/title pages now double-spaced
%
% v1.64 Mar 20 2016
%     - First paragraphs after each sectional heading is
%       now indented

% v1.63 Mar 3 2016
%     - Moved List of Publications to after Appendices
%     - References heading has 5cm top margin

%  v1.62 Feb 17 2016
%     - Turn off page numbering for appendices in the ToC
%     - Appendix chapters are un-bold in the ToC now
%
% v1.61 Sep 22 2015
%     - Abstract titles are now 12pt (irrespective of the font being used)
%     - Default bibliography style is now APA
%       (with \RequirePackage[natbibapa]{apacite})
%
% v1.6 Aug 21 2015
%     - (Based on feedback from Puan Nurul Hasyda of IPS)
%     - Chapters are now in bold in the TOC.
%     - Eliminated vertical spaces at the top and bottom of cover and
%       title page.
%     - Apparently if using Times, must be in 12pt. So no more 11pt,
%       12pt choices. But sans serif (Arial) is allowed again, but must
%       be at 11pt.
%     - Instead, the "times" option loads mathptmx 12pt as body font,
%       and "arial" option loads "helvet" 11pt as body font.
%
% v1.5 Apr 16 2015
%     - Thesis titles in abstracts are now double-spaced if >1 lines, and
%       are smaller than the "Abstract"/"Abstrak" keywords.
%
% v1.4 Apr 14 2015
%     - Now works with updated acronym package.
%     - Finally using bahasam babel option (and removed the french language,
%       couldn't remember why I didn't remove it after some personal tests).
%     - Page numbers in content lists are now unbold (requested by IPS).
%
% v1.3 May 12 2012
%     - Refined spacing in bibliography list and after chapter heading
%     - Fixed font size in cover (from 17pt to 18pt)
%     - User can now specify preferred document body font size (11pt or 12pt)
%        as a class option. The default is 11pt.
%     - Changed font sizes of section, subsection headings to \normalsize.
%
%  v1.2 Nov 15 2007
%     - New thesis format in force as of 16 November 2007!
%     - There is now no times, helvet or cm options. Body font is now
%       Times Roman for serif type, Helvetica for sans serif type, and
%       TX Typewriter for typewriter font.
%     - Bottom page margin increased to 4.0cm.
%     - Page number now placed 1.0cm from botom edge of paper.
%     - All text on title page are now in bold.
%     - Title of bibliography/references is now "REFERENCES".
%     - Introduced new "plate" float and "List of Plates".
%     - List of Publications is now moved to *after* appendices.
%     - Re-formatting of the table of contents:
%       * all text are unformatted (chapter headings are neither in bold nor
%         using larger font size)
%       * dotted leader between titles and page number
%       * larger separation between (sub)section numbers and titles to resemble
%         the sample in the guidelines more closely
%     - Fixed bug in definition of \chapter to honour short title in
%       \chapter[short title]{full title}
%
%  v1.1 May 20 2007
%     - re-definition of \chapter, therefore eliminating the need for \mychapter.
%       Various GUI LaTeX editors will now recognise the chapter divisions.
%     - added \RequirePackage[T1]{fontenc}.
%     - updated quickstart.tex: removed note about \mychapter.
%
%  v1.0 June 1 2006
%     - usmthesis is now provided as a document class (usmthesis.cls).
%     - introduced class options for:
%         * body font scheme (on suggestion by Dr Dhanesh)
%         * style of URLs
%         * whether figures are framed by default (on feedback from Anusha and
%           Seng Soon).
%     - updated usmthesis.tex to describe how to do common tasks in LaTeX
%     - new quickstart.tex to describe how to use usmthesis.cls and accompanying
%       template files.
%
%  v0.1 Nov 29 2005
%     - first stab at usmthesis as a package (usmthesis.sty).
%
%
% Available class options
% (e.g., \documentclass[urlsf,boxfig,txtt]{usmthesis}
%
%             *** choose only one from each category ***
%
% urlrm, urlsf, urltt
%   determines the font family to be used for typesetting URLs -- You may
%   want to select one that shows up URLs clearly with respect to the main
%   thesis font that you chose.  urlrm will cause URLs to be typeset with a
%   serif font (CM Roman or Times Roman); urlsf with a sans serif
%   font (CM Serif or Helvetica); and urltt with a typewriter font (CM
%   Typewriter or Courier).
%   The default is urltt.
%
% boxfig, noboxfig
%   determines if figures will have a visible bounding box.
%
%   The default is noboxfig.

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{usmthesis}[2016/03/22 version 1.65]

\RequirePackage{silence}
\WarningFilter{float}{Can't redefine}

% switches for URL font
\newif\if@urltt
\newif\if@urlsf
\newif\if@urlrm

% switch for border around figures
\newif\if@boxfig

% Added August 2015: times or arial as thesisfont
\newif\if@timesthesis
\DeclareOption{times}{\@timesthesistrue}
\DeclareOption{arial}{\@timesthesisfalse}

\DeclareOption{urltt}{\@urltttrue\@urlrmfalse\@urlsffalse}
\DeclareOption{urlsf}{\@urlsftrue\@urlttfalse\@urlrmfalse}
\DeclareOption{urlrm}{\@urlrmtrue\@urlttfalse\@urlsffalse}

\DeclareOption{boxfig}{\@boxfigtrue}
\DeclareOption{noboxfig}{\@boxfigfalse}

%\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\DeclareOption{10pt}{\OptionNotUsed}
\DeclareOption{11pt}{\OptionNotUsed}
\DeclareOption{12pt}{\OptionNotUsed}
%\DeclareOption{11pt}{\PassOptionsToClass{\CurrentOption}{book}}
%\DeclareOption{12pt}{\PassOptionsToClass{\CurrentOption}{book}}

\ExecuteOptions{times,urltt,noboxfig}
\ProcessOptions

\RequirePackage{etoolbox}
\RequirePackage[utf8]{inputenc}
\RequirePackage[T1]{fontenc}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sets up the fonts: Load the Times, Helvetica and
%% Courier (or TX Typewriter) fonts for serif, sans serif
%% and typewriter respectively.
%% (Helvetica is loaded in lieu of Arial, which is
%% commercial.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{mathptmx}
\RequirePackage{helvet}
\renewcommand{\ttdefault}{lmtt}

\if@timesthesis
  \LoadClass[a4paper,oneside,12pt]{book}
\else
  \LoadClass[a4paper,oneside,11pt]{book}
  \renewcommand\familydefault{\sfdefault}
\fi

\RequirePackage{ifpdf}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Page Boundaries and graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{graphicx}
\RequirePackage[a4paper,top=25mm,bottom=25mm,left=40mm,right=25mm,footskip=15mm]{geometry}

\ifpdf\pdfcompresslevel=9\else\fi

\RequirePackage[listofformat=parens]{subfig}
%\RequirePackage{subfigure}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% When "bahasa" is specified last, all the sectioning
% commands become Malay e.g. "Bab",
% "Dafter Isi" etc.  "english" is the main language
% here. Babel now supports bahasam.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage[bahasam,english]{babel}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default spacing is double.  In your actual content,
% whenever you need a paragraph to be single-spaced,
% do this:
%
% \begin{singlespace}
% This is your paragraph using single-spacing.
% \end{singlespace}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{setspace}
\linespread{2}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These packages will help you in typesetting
% mathematical symbols, formulae and equations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{amsmath,amssymb,textcomp}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For including the TOC in the TOC.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{tocbibind}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formatting of the various sectioning headings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{titlesec}
% chapters
%% Aug 2015: Chapter heading sizes must be FIXED at 14pt, according to IPS.
\titleformat{\chapter}[display]{\filcenter\bfseries}{\fontsize{14pt}{16pt}\selectfont\MakeUppercase{\chaptername{} \thechapter}}{-1.7pc}{\begin{singlespace}\bfseries\fontsize{14pt}{16pt}\selectfont\MakeUppercase}[\end{singlespace}]

%\titlespacing*{\chapter}{0pt}{*-4.5}{*6}
\titlespacing*{\chapter}{0pt}{*-4.5}{*4}
% section e.g. 4.1
\titleformat{\section}[hang]{\bfseries\normalsize}{\thesection}{6pt}{}
\titlespacing*{\section}{0pt}{*3}{*-3}
% subsection e.g. 4.1.1
\titleformat{\subsection}[hang]{\bfseries\normalsize}{\thesubsection}{6pt}{}
\titlespacing*{\subsection}{0pt}{*3}{*-3}
% subsubsection e.g. 4.1.1(a)
\renewcommand\thesubsubsection{\thesubsection(\alph{subsubsection})}
\titleformat{\subsubsection}[hang]{\bfseries}{\thesubsubsection}{6pt}{}
\titlespacing*{\subsubsection}{0pt}{*3}{*-3}
%% 20 Mar 2016: first paragraphs need to be indented
\RequirePackage{indentfirst}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formatting of the TOC, LOF and LOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage[titles,subfigure]{tocloft}
\setcounter{lofdepth}{2}
\setcounter{secnumdepth}{3}
\setcounter{tocdepth}{3}

% Introducing: List of Plates
\newcommand{\listplatename}{List of Plates}
\newlistof[chapter]{plate}{lop}{\listplatename}
\newcommand\listofplates{\begingroup
\tocfile{\listplatename}{lop}
\endgroup}

\renewcommand\cftdotsep{1}
\cftsetrmarg{2.6pc plus1fil}
\cftsetpnumwidth{2pc}

\setlength\cftbeforepartskip{3pc}
\setlength\cftbeforechapskip{1pc}
\setlength\cftbeforesecskip{1pc}
\setlength\cftbeforesubsecskip{1pc}
\setlength\cftbeforesubsubsecskip{1pc}
\setlength\cftbeforefigskip{1pc}
\setlength\cftbeforesubfigskip{1pc}
\setlength\cftbeforetabskip{1pc}
\setlength\cftbeforeplateskip{1pc}

\cftsetindents{sec}{0pt}{2.5pc}
\cftsetindents{subsec}{2.5pc}{3.05pc}
\cftsetindents{subsubsec}{5.45pc}{4.05pc}

\renewcommand\cftpartfont{}
\renewcommand\cftpartpagefont{\normalfont\normalsize}
\renewcommand\cftpartdotsep{1}
\renewcommand\cftpartleader{\normalfont\cftdotfill{\cftpartdotsep}}

\renewcommand\cftchapfont{}
\renewcommand\cftchappagefont{\normalfont\normalsize}
\renewcommand\cftchapdotsep{1}
\renewcommand\cftchapleader{\normalfont\cftdotfill{\cftchapdotsep}}
\renewcommand\cftchappresnum{\MakeUppercase{\chaptername} }
\renewcommand\cftchapaftersnum{\; --\; }
\newlength\tmplength

\renewcommand\cftfigpresnum{Figure }
\cftsetindents{fig}{0pt}{6.2pc}
\renewcommand\cftfigdotsep{\cftnodots}
%\renewcommand\cftfigpagefont{\bfseries}

\renewcommand\cftsubfigpresnum{Figure }
\cftsetindents{subfig}{0pt}{6.2pc}
\renewcommand\cftsubfigleader{\hfill}
%\renewcommand\cftsubfigpagefont{\bfseries}
%\renewcommand\cftfigfont{\mdseries\justifying}

\renewcommand\cfttabpresnum{Table }
\cftsetindents{tab}{0pt}{6.2pc}
\renewcommand\cfttabdotsep{\cftnodots}
%\renewcommand\cfttabpagefont{\bfseries}

\renewcommand\cftplatepresnum{Plate }
\cftsetindents{plate}{0pt}{6.2pc}
\renewcommand\cftplatedotsep{\cftnodots}
%\renewcommand\cftplatepagefont{\bfseries}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Header and Footer formatting: empty header, and
% page number at the bottom center.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\pagestyle{plain}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For customising floats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{rotfloat}
\floatstyle{plaintop}\restylefloat{table}
\if@boxfig\floatstyle{boxed}\else\floatstyle{plain}\fi
\restylefloat{figure}

\newfloat{plate}{htbp}{lop}[chapter]
\floatname{plate}{Plate}
\floatstyle{plain}\restylefloat{plate}
\def\theHplate{\arabic{plate}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For more compact-looking lists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%\RequirePackage{paralist}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package will help you in typesetting urls!
% If you just type a URL like
%    http://cs.usm.my/~abc/cgi/index.php?q=bye%20%bye
% something is going to get VERY wrong sooner or later.
% Typeset your URL like this instead:
% \url{http://cs.usm.my/~abc/cgi/index.php?q=bye%20%bye}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{url}
\if@urlrm\urlstyle{rm}\else\if@urltt\urlstyle{tt}\else\urlstyle{sf}\fi\fi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For typesetting the list of abbreviations and maybe help
% control their expansion in the text.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{acronym}
\renewcommand*{\acfsfont}[1]{\rmfamily #1}
\renewcommand*{\acsfont}[1]{\rmfamily #1}
\renewcommand*{\acffont}[1]{\rmfamily #1}

%% Apr 14, 2015: cater for new acronym.sty
%% Sep 10, 2015: make everything backward compatible
\@ifpackagelater{acronym}{2014/08/01}
{%
  % New acronym.sty works with the following
  \renewcommand*{\aclabelfont}[1]{\textbf{\acsfont{#1}}}
}
{%
  \renewenvironment{AC@deflist}[1]%
    {\ifAC@nolist%
     \else%
        \raggedright\begin{list}{}%
            {\settowidth{\labelwidth}{\textbf{\textrm{#1}}}%
            \setlength{\leftmargin}{\labelwidth}%
            \addtolength{\leftmargin}{\labelsep}%
            \renewcommand{\makelabel}{\bflabel}}%
      \fi}%
    {\ifAC@nolist%
     \else%
        \end{list}%
     \fi}
}%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table-related packages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{tabularx,colortbl,dcolumn}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 22 Sep 2015: Loads apacite by default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \RequirePackage[numbers]{natbib}
% \RequirePackage{natbib}
\RequirePackage[natbibapa]{apacite}
\addto\captionsenglish{%
\renewcommand{\bibname}{REFERENCES}
}
\setlength\bibsep{1em}
\setlength\bibhang{1.5em}
\renewcommand\bibfont{\singlespace}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiple bibliographics, used to hack up a list of
% own publications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage[resetlabels]{multibib}

%% Title of the thesis in Malay
\def\titlems#1{\gdef\@titlems{#1}}
\def\@titlems{\@latex@warning@no@line{No \noexpand\titlems given}}

%% Degree type: MSc or PhD
\def\degreetype#1{\gdef\@degreetype{#1}}
\def\@degreetype{\@latex@warning@no@line{No \noexpand\degreetype given}}

%% Month and year in which thesis was submitted
\def\submitmonth#1{\gdef\@submitmonth{#1}}
\def\@submitmonth{\@latex@warning@no@line{No \noexpand\submitmonth given}}
\def\submityear#1{\gdef\@submityear{#1}}
\def\@submityear{\@latex@warning@no@line{No \noexpand\submityear given}}

%% overrides \chapter, to allow for all-caps in text but as-is in TOC
%\newcommand\mychapter[1]{\chapter[\texorpdfstring{\uppercase{#1}}{#1}]{#1}}

% May 20: re-defined \@chapter so \mychapter is not needed anymore. This
% allows IDEs such as TeXnicCenter to properly recognise the chapters and
% display them correctly in the Structure/Navigator.
\def\@chapter[#1]#2{\ifnum \c@secnumdepth >\m@ne
                       \if@mainmatter
                         \refstepcounter{chapter}%
                         \typeout{\@chapapp\space\thechapter.}%
                         \addcontentsline{toc}{chapter}%
                                   {\protect\numberline{\thechapter}\texorpdfstring{\uppercase{#1}}{#1}}%
                       \else
                         \addcontentsline{toc}{chapter}{#1}%
                       \fi
                    \else
                      \addcontentsline{toc}{chapter}{\texorpdfstring{\uppercase{#1}}{#1}}%
                    \fi
                    \chaptermark{\texorpdfstring{\uppercase{#1}}{#1}}%
                    \addtocontents{lof}{\protect\addvspace{10\p@}}%
                    \addtocontents{lot}{\protect\addvspace{10\p@}}%
                    \addtocontents{lop}{\protect\addvspace{10\p@}}%
                    \if@twocolumn
                      \@topnewpage[\@makechapterhead{\texorpdfstring{\uppercase{#2}}{#1}}]%
                    \else
                      \@makechapterhead{\texorpdfstring{\uppercase{#2}}{#1}}%
                      \@afterheading
                    \fi
}



%% Creates the first two cover pages
\newcommand\makecover{%
\begin{titlepage} \centering\fontsize{18pt}{20pt}\selectfont\bfseries{%
%\vspace*{10mm}
%
{\fontsize{18pt}{20pt}\selectfont\MakeUppercase{\@title}\par}
%
\vspace*{\stretch{1}}

\MakeUppercase{\@author}\par

\vspace*{\stretch{1}}

Certificate in Quantitative Finance\par

\@submityear
%
%\vspace*{20mm}
}

\end{titlepage}


\thispagestyle{empty}{\centering\bfseries
%\vspace*{15mm}

{\fontsize{18pt}{20pt}\selectfont\MakeUppercase{\@title}\par}

\vspace*{\stretch{2}}

{\Large by}\par

\vspace*{\stretch{2}}

\MakeUppercase{\fontsize{18pt}{20pt}\selectfont\@author}\par

\vspace*{\stretch{2}}

\begin{singlespace}
\fontsize{16pt}{18pt}\selectfont Final report submitted in fulfilment of the requirements\\
for \fontsize{16pt}{18pt}\selectfont \@degreetype\par
\end{singlespace}

\vspace*{\stretch{2}}

{\fontsize{18pt}{20pt}\selectfont\@submitmonth \ \@submityear }\par

%\vspace*{20mm}
}
\clearpage}


%% English environment
\newenvironment{EnAbstract}
  { \titleformat{\chapter}[display]{\filcenter\bfseries}{\fontsize{12pt}{14pt}\selectfont\MakeUppercase{\chaptername{} \thechapter}}{-1.7pc}{\linespread{2}\bfseries\fontsize{12pt}{14pt}\selectfont\MakeUppercase}
    \chapter[Abstract]{\@title}%
    \begin{center}{\fontsize{12pt}{14pt}\bfseries\selectfont ABSTRACT}\end{center}\par\indent}
  {\ignorespacesafterend}

%% Malay abstract environment
\newenvironment{MsAbstract}
  {\begin{otherlanguage}{bahasam}%
    \titleformat{\chapter}[display]{\filcenter\bfseries}{\fontsize{12pt}{14pt}\selectfont\MakeUppercase{\chaptername{} \thechapter}}{-1.7pc}{\linespread{2}\bfseries\fontsize{12pt}{14pt}\selectfont\MakeUppercase}
    \chapter[Abstrak]{\@titlems}%
    \begin{center}{\fontsize{12pt}{14pt}\bfseries\selectfont ABSTRAK}\end{center}\par\indent}
  {\end{otherlanguage}\ignorespacesafterend}

\apptocmd{\frontmatter}{%
  \addtocontents{toc}{\protect\raggedright\protect\sloppy}
  \addtocontents{toc}{\protect\hfill{\protect\bfseries Page}\\}
  \addtocontents{lof}{\protect\raggedright\protect\sloppy}
  \addtocontents{lof}{\protect\hfill{\protect\bfseries Page}\\}
  \addtocontents{lot}{\protect\raggedright\protect\sloppy}
  \addtocontents{lot}{\protect\hfill{\protect\bfseries Page}\\}
  \addtocontents{lop}{\protect\raggedright\protect\sloppy}
  \addtocontents{lop}{\protect\hfill{\protect\bfseries Page}\\}
  \renewcommand\contentsname{Table of Contents}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Inserts a "Front Matter" bookmark if using pdfLaTeX
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  \ifpdf
      \pdfbookmark[-1]{Front Matter}{front}
  \else\fi
  \renewcommand{\thechapter}{\arabic{chapter}}
}{}{}

\pretocmd{\mainmatter}{%
  \addtocontents{toc}{\protect\setlength{\protect\cftbeforechapskip}{2.5pc}}
  \addtocontents{toc}{\protect\cftpagenumbersoff{chap}}
  \addtocontents{toc}{\protect\renewcommand\protect\cftchapfont{\protect\bfseries}}
}{}{}

\apptocmd{\mainmatter}{%
  \ifpdf\pdfbookmark[-1]{Main Chapters}{main}\else\fi
}{}{}

\pretocmd{\bibliography}{%
  \addtocontents{toc}{\protect\cftpagenumberson{chap}}
  \begin{singlespace}
}{}{}

\apptocmd{\bibliography}{%
  \end{singlespace}
  \clearpage
}{}{}

\pretocmd{\appendix}{%
  \clearpage
  \addtocontents{toc}{\protect\cftpagenumbersoff{part}}
  \addtocontents{toc}{\protect\cftpagenumbersoff{chapter}}
  \addtocontents{toc}{\protect\setlength{\protect\cftbeforechapskip}{1pc}}
}{}{}

\apptocmd{\appendix}{
  \addtocontents{toc}{\protect\renewcommand\protect\cftchapfont{\protect\mdseries}}
}
{}{}

\AtBeginDocument{
\settowidth{\tmplength}{\bfseries\MakeUppercase{\chaptername}\space9\cftchapaftersnum}
\setlength\cftchapnumwidth{\tmplength}
\apptocmd{\appendix}{%
  \phantomsection\addcontentsline{toc}{part}{\texorpdfstring{\bfseries\uppercase{Appendices}}{Appendices}}
  \thispagestyle{empty}
  \vspace*{\stretch{1}}
    \begin{center}
      {\Huge\bfseries APPENDICES}
    \end{center}
  \vspace*{\stretch{2}}
  \addtocontents{toc}{\protect\renewcommand{\protect\chaptername}{\appendixname}}
  \settowidth{\tmplength}{\bfseries\MakeUppercase{\appendixname}\space N\cftchapaftersnum}
  \addtocontents{toc}{\protect\setlength{\cftchapnumwidth}{\the\tmplength}}
  \renewcommand\chaptername{\appendixname}
}{}{}
}
