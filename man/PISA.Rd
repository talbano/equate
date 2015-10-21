\name{PISA}
\alias{PISA}
\docType{data}
\title{Programme for International Student Assessment 2009 USA Data}
\description{
This dataset contains scored cognitive item response data from the 2009 administration of the Programme for International Student Assessment (PISA), an international study education systems. The data, along with the license under which they are released, are available online at \url{http://www.oecd.org/pisa/}.
}
\usage{PISA}
\format{
\code{PISA} is a \code{list} containing four elements. The first, \code{PISA$students}, is a \code{data.frame} containing 233 variables across 5233 individuals, with one row per individual. All but one variable come from the USA PISA data file "INT_COG09_S_DEC11.txt". The remaining variable, language spoken at home, has been merged in from the student questionnaire file "INT_STQ09_DEC11.txt". Variable names match those found in the original files:
\describe{
  \item{\code{stidstd}}{
  	Unique student ID (one for each of the 5233 cases);
  }
  \item{\code{schoolid}}{
  	School ID (there are 165 different schools);
  }
  \item{\code{bookid}}{
  	ID for the test booklet given to a particular student, of which there were 13;
  }
  \item{\code{langn}}{
  	Student-reported language spoken at home, with 4466 students reporting English (indicated by code 313), 484 students reporting Spanish (with code 156) and 185 students reporting "another language" (code 859);
  }
  \item{\code{m033q01} to \code{s527q04t}}{
  	Scored item-response data across the 189 items included in the general cognitive assessment, described below; and
  }
  \item{\code{pv1math} to \code{pv5read5}}{
  	PISA scale scores, referred to in the PISA technical documentation as "plausible values".
  }
}
Next, \code{PISA$booklets} is a \code{data.frame} containing 4 columns and 756 rows and describes the 13 general cognitive assessment booklets. Variables include:
\describe{
	\item{\code{bookid}}{
		The test booklet ID, as in \code{PISA$students};
	}
	\item{\code{clusterid}}{
		ID for the cluster or item subset in which an item was placed; items were fully nested within clusters; however, each item cluster appeared in four different test booklets;
	}
	\item{\code{itemid}}{
		Item ID, matching the columns of \code{PISA$students}; each item appears in \code{PISA$booklets} four times, once for each booklet; and
	}
	\item{\code{order}}{
		The order in which the cluster was presented within a given booklet.
	}
}
\code{PISA$items} is a \code{data.frame} containing 4 columns and 189 rows, with one row per item. Variables include:
\describe{
	\item{\code{itemid}}{
		Item ID, as in \code{PISA$booklets}
	}
	\item{\code{clusterid}}{
		Cluster ID, as in \code{PISA$booklets}
	}
	\item{\code{max}}{
		Maximum possible score value, either 1 or 2 points, with dichotomous scoring (max of 1) used for the majority of items; and
	}
	\item{\code{subject}}{
		The subject of an item, equivalent to the first character in \code{itemid} and \code{clusterid}.
	}
}
Finally, \code{PISA$totals} is a list of 13 \code{data.frame}s, one per booklet, where the columns correspond to total scores for all students on each cluster for the corresponding booklet. These total scores were calculated using \code{PISA$students} and \code{PISA$booklets}. Elements within the \code{PISA$totals} list are named by booklet, and the columns in the \code{data.frame} are named by cluster. For example, \code{PISA$totals$b1$m1} contains the total scores on cluster M1 for students taking booklet 1.
}
\source{
OECD (2012). PISA 2009 Technical Report, PISA, OECD Publishing. \url{http://dx.doi.org/10.1787/9789264167872-en}

Addition information can be found at the PISA website: \url{http://www.oecd.org/pisa/}
}
\keyword{datasets}
