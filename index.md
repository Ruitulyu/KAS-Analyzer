## KAS-pipe2: a user-friendly toolkit for exploring KAS-seq data

### Author: Ruitu Lyu, Chemistry department, the University of Chicago

----------------------------------------

<img src="https://github.com/Ruitulyu/KAS-pipe2/blob/main/image/KAS-pipe2.jpg"  />

### Introduction

KAS-pipe2 is a collection of command line tools specifically developped for exploring KAS-seq or strand-specific (sp)KAS-seq data, including the basic processing tools for quality control, reference genome index, raw reads mapping, and heatmaps and summary plots. KAS-pipe2 also includes many novel features and completely new frame design compared to KAS-pipe. e.g. time-courese(TC) KAS-seq differential analysis, R-loop identification (only for spKAS-seq), ss-enhancers, motif, index calculation, termination length and so on.

KAS-pipe2 is still on active development and the source code is hosted on [GitHub](https://github.com/Ruitulyu/KAS-pipe2).

### Installation

**Install by cloning KAS-pipe2 git repository on github:**

You can install KAS-pipe2 on command line (linux/mac) by cloning git repository on github:

	$ git clone https://github.com/Ruitulyu/KAS-pipe2.git
	$ cd KAS-pipe2
	$ bash ./setup.sh
	
	# If anaconda or miniconda was not installed on your system.
	$ KAS-pipe2 install -conda
	
	# Install conda 'KAS-pipe2' environment. 
	$ KAS-pipe2 install -KAS-pipe2
	
	# Activate conda 'KAS-pipe2' environment.
	$ KAS-pipe2 activate
	

You can use the [editor on GitHub](https://github.com/Ruitulyu/KAS-pipe2/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.
<script type="text/javascript" src="//rf.revolvermaps.com/0/0/6.js?i=5f2zpl0mkwd&amp;m=7&amp;c=e63100&amp;cr1=ffffff&amp;f=arial&amp;l=0&amp;bv=90&amp;lx=-420&amp;ly=420&amp;hi=20&amp;he=7&amp;hc=a8ddff&amp;rs=80" async="async"></script>

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/Ruitulyu/KAS-pipe2/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
