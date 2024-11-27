import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { Scatter } from 'react-chartjs-2';
import Papa from 'papaparse';

const DataGenerator = () => {
  const [inputValue, setInputValue] = useState('');
  const [labels, setLabels] = useState([]);
  const [selectedLabel, setSelectedLabel] = useState('');
  const [chartData, setChartData] = useState({});
  const [error, setError] = useState('');

  const yLabels = [
    'HXD13', 'HXD12', 'HXD11', 'HXD10', 'HXD9', 'HXD8', 'HXD4', 'HXD3', 'HXD1', 
    'HXC13', 'HXC12', 'HXC11', 'HXC10', 'HXC9', 'HXC8', 'HXC6', 'HXC5', 'HXC4', 
    'HXB13', 'HXB9', 'HXB8', 'HXB7', 'HXB6', 'HXB5', 'HXB4', 'HXB3', 'HXB2', 
    'HXB1', 'HXA13', 'HXA11', 'HXA10', 'HXA9', 'HXA7', 'HXA6', 'HXA5', 'HXA4', 
    'HXA3', 'HXA2', 'HXA1'
  ];

  const aminoAcidLength = [
    [39, 335], [38, 376], [37, 443], [36, 320], [35, 270], [34, 233], [33, 230], [32, 272], 
    [31, 410], [30, 313], [29, 388], [28, 301], [27, 356], [26, 431], [25, 251], 
    [24, 269], [23, 224], [22, 217], [21, 243], [20, 250], [19, 284], [18, 264], 
    [17, 222], [16, 235], [15, 242], [14, 260], [13, 342], [12, 304], [11, 282], 
    [10, 330], [9, 328], [8, 432], [7, 255], [6, 290], [5, 352], [4, 340], 
    [3, 338], [2, 270], [1, 343]
  ];

  useEffect(() => {
    fetchLabels();
  }, []);

  const fetchLabels = async () => {
    try {
      const response = await axios.get('http://localhost:5000/labels');
      setLabels(response.data.labels);
    } catch (error) {
      console.error('Error fetching labels!', error);
      setError('Error fetching labels!');
    }
  };

  const handleGenerateData = async () => {
    try {
      const response = await axios.post('http://localhost:5000/generate-data', { input: inputValue });
      const { label } = response.data;
      setLabels([...labels, label]);
      setInputValue('');
    } catch (error) {
      console.error('There was an error generating the data!', error);
      setError('There was an error generating the data!');
    }
  };

  const handleDownloadData = async (label) => {
    try {
      const response = await axios.get(`http://localhost:5000/data/${label}.csv`, {
        responseType: 'blob',
      });
      const url = window.URL.createObjectURL(new Blob([response.data]));
      const link = document.createElement('a');
      link.href = url;
      link.setAttribute('download', `${label}.csv`);
      document.body.appendChild(link);
      link.click();
      link.parentNode.removeChild(link);
    } catch (error) {
      console.error('Error downloading data!', error);
      setError('Error downloading data!');
    }
  };

const parseCsvData = (csvData, label) => {
  const arrayColumns = ['Residues', 'Weights', 'Amino Acids'];

  const parsedCsv = Papa.parse(csvData, {
    header: true,
    skipEmptyLines: true,
  });

  const data = parsedCsv.data;

  const parseArray = (value) => {
    try {
      return JSON.parse(value.replace(/'/g, '"'));
    } catch (e) {
      return value ? [value] : []; // Return an array with the value if parsing fails
    }
  };

  const expandRow = (row, arrayColumns) => {
    const maxLen = Math.max(...arrayColumns.map(col => row[col].length));
    let expandedRows = [];

    for (let i = 0; i < maxLen; i++) {
      let newRow = {};
      for (let key in row) {
        if (arrayColumns.includes(key)) {
          newRow[key] = row[key][i] !== undefined ? row[key][i] : row[key][row[key].length - 1];
        } else {
          newRow[key] = row[key];
        }
      }
      expandedRows.push(newRow);
    }

    return expandedRows;
  };

  let processedData = [];

  data.forEach(row => {
    arrayColumns.forEach(col => {
      row[col] = parseArray(row[col]);
    });
    const expandedRows = expandRow(row, arrayColumns);
    processedData = processedData.concat(expandedRows);
  });

  const dataEntries = processedData
    .filter(row => parseFloat(row.Residues) !== 0) // Filter out rows with residue value of 0
    .map(row => ({
      x: parseFloat(row.Residues),
      y: parseFloat(row.Weights),
      identifier: row['Amino Acids'], // Include amino acid as the identifier
    }));

  console.log("Parsed Data Entries:", dataEntries);

  return {
    label: label,
    data: dataEntries,
    pointRadius: 3,
    backgroundColor: 'rgba(75, 192, 192, 0.5)', // Use rgba for alpha
    borderColor: 'rgba(75, 192, 192, 0.7)', // Border with alpha
    pointLabelFontSize: 14, // Optional: Font size for point labels
  };
};

const handlePlotData = async (label) => {
  try {
    const response = await axios.get(`http://localhost:5000/data/${label}.csv`);
    const { data } = response;
    const parsedData = parseCsvData(data, label);

    if (parsedData) {
      setSelectedLabel(label);
      const horizontalLineDatasets = generateHorizontalLines(aminoAcidLength);
      setChartData({
        datasets: [parsedData, ...horizontalLineDatasets],
      });
      logScatterPlotData(parsedData.data);
    } else {
      console.error('Data for the selected label not found!');
      setError('Data for the selected label not found!');
    }
  } catch (error) {
    console.error('Error plotting data!', error);
    setError('Error plotting data!');
  }
};

  

  const generateHorizontalLines = (linesData) => {
    return linesData.map(([y, xLimit], index) => ({
      label: `Line ${index + 1}`,
      data: [{ x: 0, y }, { x: xLimit, y }],
      backgroundColor: 'rgba(0, 0, 0, 0)',
      borderColor: 'rgba(0, 0, 0, 0.5)',
      borderWidth: 1,
      pointRadius: 0,
      showLine: true,
      fill: false,
      skipLegend: true,
    }));
  };

  const logScatterPlotData = (data) => {
    console.log("Scatter Plot Data:");
    data.forEach(point => {
      console.log(`Protein: ${point.protein_name}, Residue: ${point.x}, Ranked Weight: ${point.y}`);
    });
  };

  const options = {
    animation: false,
    scales: {
      x: {
        type: 'linear',
        position: 'bottom',
        title: {
          display: true,
          text: 'Amino Acid Position',
        },
        ticks: {
          padding: 10,
        }
      },
      y: {
        title: {
          display: true,
          text: 'Protein',
        },
        ticks: {
          min: 1,
          max: yLabels.length,
          stepSize: 1,
          autoSkip: false,
          padding: 10,
          callback: function(value) {
            return yLabels[value - 1];
          }
        }
      }
    },
    layout: {
      padding: {
        left: 10,
        right: 10,
        top: 10,
        bottom: 10
      }
    },    
    plugins: {
      tooltip: {
        callbacks: {
          label: function(tooltipItem) {
            const dataset = chartData.datasets[tooltipItem.datasetIndex];
            const dataPoint = dataset.data[tooltipItem.dataIndex];
            return `${dataPoint.protein_name}: (x: ${dataPoint.x}, y: ${dataPoint.y})`;
          }
        }
      },
      legend: {
        display: true,
        onClick: () => null,
        labels: {
          filter: (legendItem) => {
            return !legendItem.text.startsWith('Line');
          }
        }
      },
    }
  };

  return (
    <div>
      <h2>Data Generator</h2>
      {error && <p style={{ color: 'red' }}>{error}</p>}
      <div>
        <input
          type="text"
          value={inputValue}
          onChange={(e) => setInputValue(e.target.value)}
          placeholder="Enter input value"
        />
        <button onClick={handleGenerateData}>Generate Data</button>
      </div>
      <div>
        <h3>Labels:</h3>
        <ul>
          {labels.map((label, index) => (
            <li key={index}>
              {label}
              <button onClick={() => handleDownloadData(label)}>Download Data</button>
              <button onClick={() => handlePlotData(label)}>Plot Data</button>
            </li>
          ))}
        </ul>
      </div>
      {selectedLabel && (
        <div style={{ marginTop: '20px' }}>
          <h3>Plot for {selectedLabel}</h3>
          <Scatter 
            data={chartData} 
            options={options}
          />
        </div>
      )}
    </div>
  );
};

export default DataGenerator;
