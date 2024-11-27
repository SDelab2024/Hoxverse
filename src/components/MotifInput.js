import React, { useState } from 'react';
import axios from 'axios';

const MotifInput = ({ setInitialDatasetsType3 }) => {
  const [motif, setMotif] = useState('');
  const [error, setError] = useState('');

  const handleSubmit = async (e) => {
    e.preventDefault();
    setError('');  // Clear previous errors

    try {
      const response = await axios.post('http://127.0.0.1:5000/api/datasetsType3', { motif });
      if (response.data.error) {
        setError(response.data.error);  // Display error message
      } else {
        setInitialDatasetsType3(response.data);
      }
    } catch (error) {
      console.error('Error fetching data:', error);
    }
  };

  return (
    <form onSubmit={handleSubmit}>
      <label>
        Enter Motif:
        <input
          type="text"
          value={motif}
          onChange={(e) => setMotif(e.target.value.toUpperCase())}
          pattern="[A-Z]+"
          title="Motif must contain only alphabetical characters"
          required
        />
      </label>
      <button type="submit">Submit</button>
      {error && <p style={{ color: 'red' }}>{error}</p>}  {/* Display error message */}
    </form>
  );
};

export default MotifInput;
