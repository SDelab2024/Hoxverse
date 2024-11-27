import React from 'react';

const LabelSelector = ({ onLabelSelected, datasets }) => {
  return (
    <div>
      <h2>Select Datasets</h2>
      <ul>
        {datasets.map((dataset) => (
          <li key={dataset.label}>
            <input
              type="checkbox"
              id={dataset.label}
              onChange={() => onLabelSelected(dataset.label)}
            />
            <label htmlFor={dataset.label}>{dataset.label}</label>
          </li>
        ))}
      </ul>
    </div>
  );
};

export default LabelSelector;
