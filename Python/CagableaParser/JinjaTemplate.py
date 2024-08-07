
# https://realpython.com/primer-on-jinja-templating/

template = """
    <html>
    <head>
    <style>
    .dotForTables {
            border:1px solid black;
            height:10px;
            width:10px;
            font-size:13px;
            border-radius: 50%;
            display:inline-block;
            }
    body {
          background-color: white;
        }

h2 {
 color: darkblue;
 text-align: left;
}
table {
  border-collapse: separate;
  border-spacing: 10px;
  /* Apply cell spacing */
}
 td {
  padding: 15px;
  text-align: center;
}


.styled-table {
    border-collapse: collapse;
    margin: 10px 0;
    font-size: 0.9em;
    font-family: sans-serif;
    min-width: 400px;
    box-shadow: 0 0 20px rgba(0, 0, 0, 0.15);
}

.styled-table thead tr {
    background-color: #009879;
    color: #ffffff;
    text-align: left;
}

.styled-table th,
.styled-table td {
    padding: 12px 15px;
}

.styled-table tbody tr {
    border-bottom: 1px solid #dddddd;
}

.styled-table tbody tr:nth-of-type(even) {
    background-color: #f3f3f3;
}

.styled-table tbody tr:last-of-type {
    border-bottom: 2px solid #009879;
}

.styled-table tbody tr.active-row {
    font-weight: bold;
    color: #009879;
}
        h1 {
          color: maroon;
          margin-left: 40px;
        }

.maingrid {
        display: grid;
    align-content: start;
        grid-template-columns: repeat(2, [col] auto ) ;
        grid-gap: 10px;
        width: 100%;
    }

   .divTable{
	display: table;
	width: 100%;

}
.divTableRow {
	display: table-row;
}
.divTableHeading {
	background-color: #EEE;
	display: table-header-group;
}
.divTableCell, .divTableHead {
    text-align: center;
	border: 1px solid #999999;
	display: table-cell;
	padding: 2px 3px;
}
.divTableHeading {
	background-color: #EEE;
	display: table-header-group;
	font-weight: bold;
}
.divTableFoot {
	background-color: #EEE;
	display: table-footer-group;
	font-weight: bold;
}
.divTableBody {
	display: table-row-group;
} 

</style>
</head>

    <body>
   <h2>Samples:</h2>
    <table style="width:100%">

    {% for onelist in samples %}
    <tr>
    {% for oneitem in onelist %}
    <td>{{oneitem}}</td> 
    {% endfor %}
    </tr>
    {% endfor %}
    </h3>
    </table>
    <div>
    <h2>Parameters</h2>
    <table style="width:50%">
    <tr><td>Min depth:</td><td>{{minDepth}}</td></tr>
    <tr><td>Min frequency:</td><td>{{minFreq}}</td></tr>
    <tr><td>Min frequency Heatmap:</td><td>{{minFreqForHeatMaps}}</td></tr>

    <tr><td>Min indel frequency:</td><td>{{minIndelFreq}}</td></tr>
    <tr><td>Range:</td><td>{{plotrg}}</td></tr>
    </table>
    </div>
    {% if heatmap is not none %}
        <h2>Heatmap</h2>
        {{ heatmap }}
    {%endif%}
    {% if clustermap is not none %}
        <div>
        <h2>Clustered Heatmap</h2>
        {{ clustermap }}
        </div>
    
    {%endif%}
    {% if variantsHeatmap is not none %}
        <div>
        <h2>Variants Heatmap</h2>
        {{ variantsHeatmap }}
        </div>
    
    
    {%endif%}
    {% if depth is not none %}
    <h2>Sequencing depths</h2>
    <h4>Coverage was calculated relative to the region of the genome theoretically covered by amplicons: ({{settings.maxAmplifiedRange[0]}} - {{settings.maxAmplifiedRange[1]}})</h4>
    {{ depth }}
    {%endif%}
    {% if violins is not none %}
    <h2>Sequencing depths frequencies</h2>
    {{ violins }}
    {%endif%}
    {% if shannon is not none %}
    <h2>Shannon Entropy Distribution</h2>
    {{ shannon }}
    {%endif%}
     {% if shannonBoxFig is not none %}
    <h2>Shannon Entropy Boxplots</h2>
    {{ shannonBoxFig }}
    {%endif%}

    {% if pepperDF is not none %} 
        <div>
        <h2>Sars-Cov-2/pepper count ratios</h2>    
        {{ pepperDF.to_html(classes='styled-table') | safe}}
        {{ pepperfig.to_html(full_html=False, include_plotlyjs='cdn') }}
        </div>
    {%endif%}
    


<h2>Variants</h2>
Variants tested: 
{% for node in nodes %}
    {{node}} / 
{% endfor %}
{% if variantData.variantPieList is not none %}
    <table>
     {% for onelist in samples %}
     <tr>
     {% for oneitem in onelist %}
     {% if variantData.variantPieList[oneitem] is not none %}
     <td>
     <b>{{oneitem}}</b><br>
        {{variantData.variantPieList[oneitem]}}
     </td>
     {% endif %}
     {% endfor %}
     </tr>
     {% endfor %}
    </table>>
{%endif%}
{{variantData.variantHistoFigHTML}}

<h1>Detailed Counts</h1>
{% set meanadjustedCauseParentFlag = '&#9757;' %}
{% set requiredMutationFlag = '&#9967;' %}

{% if settings.verbose %}
    <br>For certain variants Means are corrected for other variants that are not subvariants but rather should be. E.g. for BA.2 mean, BA.4 and BA.5 is substracted since they share almost all BA.2 mutations --> means shown do not necessarily reflect the means of the values of the individual mutations. 
    <br><br>{{meanadjustedCauseParentFlag}} before mean indicates that mean was corrected since sum of child freqs exceeded parent freq <br> 
    <br>red values indicates outliers.<br> 
    <br>{{requiredMutationFlag}} in front of a mutation means this mutation is important </br>
    If certain positions for a variant are positive but in total less than the minimal number of positions (defined for each variant) are positive, the mean is set to 0 and printed in blue.
    <br>   
{%endif%}  
    {% for oneVariant in variantData.variantsDetailedCounts %}    
    <h2>{{oneVariant}}</h2>
    {% if nodes[oneVariant] is not none %}
    If less than {{nodes[oneVariant].minmutsforpass |string}} have frequencies  > 0 mean is set to 0 (mean printed blue)
    {% endif %}
    {% if nodes[oneVariant].requiredMutations is not none %}
    {% if nodes[oneVariant].minstarredmutsforpass != 0 %}
    <br> If less than {{nodes[oneVariant].minstarredmutsforpass |string}} (marked with '+') have above 0 frequencies, mean is set to 0 (mean printed blue)</br>
    {% endif %}
    {% endif %}
    {% if nodes[oneVariant] is not none %}
    {% if nodes[oneVariant].comment is not none %}
    <br>{{nodes[oneVariant].comment}}
    {% endif %}
    {% endif %}
 
    <table class="styled-table">
    <tr>
    <th></th>
    {% for oneSample in variantData.variantDataMeans %}
    <th><b>{{oneSample}}</b> </th>
    {% endfor %}
    </tr> <!–– end divTableRow -->
    {% for onemut in variantData.variantsDetailedCounts[ oneVariant].index %}
        <tr>
        <th><b>
        {% if nodes[oneVariant].testIfMutIsMandatory(onemut)%}{{requiredMutationFlag}}{% endif %}{{onemut}}</b>
        {{variantData.variantsDetailedCounts[oneVariant].loc[onemut]['dots']}}
        </th>
        
        <!–– print counts -->
        {% for oneSample in variantData.variantDataMeans %}
            <td>
            {% if variantData.variantsdetailedCountsMask[oneVariant][oneSample].loc[onemut] == false %}
                <font color="red">{{"%.3f"|format(variantData.variantsDetailedCounts[oneVariant][oneSample].loc[onemut])}}</font>
            {% else %}
                {{"%.3f"|format(variantData.variantsDetailedCounts[oneVariant][oneSample].loc[onemut])}}
            {% endif %}
            </td>
        {% endfor %}
        </tr> 

    {% endfor %}    <!–– end loop over mutations -->
     <!–– add row with means -->
        <tr>
        <th><b>Mean</b></th>
        {% for oneSample in variantData.variantDataMeans %}
            {% if variantData.variantDataMeans[oneSample].loc[oneVariant].meanSetToZero %}
                {% set col = 'blue' %}             
            {% else %}
                {% set col = 'black' %}
            {% endif %}
            {% set s = '' %}
            {% if variantData.variantDataMeans[oneSample].loc[oneVariant].meanWasCorrectedToFitParents %}
                {% set s = meanadjustedCauseParentFlag %}
            {% endif %}
            <td><b><font color={{col}}>{{s}}{{"%.3f"|format(variantData.variantDataMeans[oneSample].loc[oneVariant].mean) }}</font></b></td> 
        {% endfor %}
        </tr>
             <!–– add row with SE -->
        <tr>
        <th><b>SE</b></th>
        {% for oneSample in variantData.variantDataMeans %}
          <td>  {{"%.3f"|format(variantData.variantDataMeans[oneSample].loc[oneVariant].se) }}</td>
        {% endfor %}
        </tr>
         <tr>
        <th></th>
        {% for oneSample in variantData.variantDataMeans %}
          <td>  {{variantData.variantDataMeans[oneSample].loc[oneVariant].box}}</td>
        {% endfor %}
        </tr>

    </table> <!–– end one variant table -->
     {% endfor %}

<br>

    </body>
    </html>
"""
