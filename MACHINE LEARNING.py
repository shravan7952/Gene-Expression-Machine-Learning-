import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import joblib
import seaborn as sns
import os
from sklearn.metrics import confusion_matrix

# List of cancer data files
cancer_files = [
    'breast cancer (BRACA).tsv', 'Adenoid cystic carcinoma (ACC).tsv', 'Bladder urothelial carcinoma (BLCA).tsv',
    'Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (CESE).tsv', 'Cholangiocarcinoma (CHOL).tsv',
    'Diffuse large B-cell lymphoma (DLBCL).tsv', 'Esophageal cancer (ESCA).tsv', 'glioblastoma multiforme (GBM).tsv',
    'Head and neck squamous cell carcinoma (HNSC).tsv', 'Kidney renal clear cell carcinoma (KIRC).tsv',
    'Kidney renal papillary cell carcinoma (KIRP).tsv', 'Acute myeloid leukemia (LAML).tsv',
    'Pediatric low-grade gliomas (LGG).tsv', 'Liver Hepatocellular Carcinoma(LIHC).tsv', 'Lung adenocarcinoma (LUAD).tsv',
    'lung squamous cell carcinoma (LUSC).tsv', 'mesothelioma cancer (MESO).tsv', 'Ovarian cancer (OV).tsv',
    'Pancreatic adenocarcinoma (PAAD).tsv', 'Pheochromocytoma and paraganglioma (PCPGs).tsv',
    'Prostate Adenocarcinoma (PRAD).tsv', 'Rectum Adenocarcinoma (READ).tsv', 'Sarcoma Cancer (SARC).tsv',
    'Skin cutaneous melanoma (SKCM).tsv', 'Stomach adenocarcinoma (STAD).tsv', 'Tenosynovial giant cell tumors (TGCTs).tsv',
    'Thyroid carcinoma (THCA).tsv', 'Thymoma and Thymic Carcinoma(THYM).tsv', 'Uterine corpus endometrial carcinoma (UCEC).tsv',
    'Uterine Carcinosarcoma (UCS).tsv', 'Uveal melanoma (UVM).tsv', 'Colon adenocarcinoma (COAD).tsv',
    'Kidney Chromophobe (KICH).tsv'
]

# Load data from files
combined_data = []
for cancer_file in cancer_files:
    df = pd.read_csv(cancer_file, sep='\t', index_col=0)
    cancer_type = cancer_file.split('.')[0]
    df['Cancer_Type'] = cancer_type
    combined_data.append(df)

# Combine all data into one DataFrame
final_data = pd.concat(combined_data, axis=0)
final_data.reset_index(drop=True, inplace=True)
final_data.fillna(0, inplace=True)

# Prepare features and labels
X = final_data.drop(columns=['Cancer_Type']).values
y = final_data['Cancer_Type'].values

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42)

# Encode labels
label_encoder = LabelEncoder()
y_train_encoded = label_encoder.fit_transform(y_train)
y_test_encoded = label_encoder.transform(y_test)

# Train the model
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train, y_train_encoded)

# Evaluate the model
y_pred_encoded = model.predict(X_test)
accuracy = accuracy_score(y_test_encoded, y_pred_encoded)
print(f"Accuracy: {accuracy * 100:.2f}%")
print(classification_report(y_test_encoded, y_pred_encoded, target_names=label_encoder.classes_))

# Feature importance
importances = model.feature_importances_
sorted_idx = importances.argsort()
plt.barh(final_data.columns[sorted_idx][-10:], importances[sorted_idx][-10:])
plt.xlabel('Feature Importance')
plt.title('Top 10 RNA genes')
plt.show()

# Save the model and label encoder
joblib.dump(model, 'cancer_classification_model.joblib')
joblib.dump(label_encoder, 'label_encoder.joblib')
print("Model and Label Encoder saved successfully!")

# Standardize the data for PCA
scaler = StandardScaler ()
X_scaled = scaler.fit_transform(X)

# Perform PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# Create a DataFrame for PCA results
pca_df = pd.DataFrame(data=X_pca, columns=['PC1', 'PC2'])
pca_df['Cancer_Type'] = y

# Plotting the PCA results
plt.figure(figsize=(10, 8))
cancer_types = label_encoder.classes_
markers = ['o', 's', 'D', '^', 'v', 'p', '*', 'H', 'X', '<', '>']

for i, cancer_type in enumerate(cancer_types):
    subset = pca_df[pca_df['Cancer_Type'] == cancer_type]
    plt.scatter(subset['PC1'], subset['PC2'], 
                label=cancer_type, 
                marker=markers[i % len(markers)], 
                edgecolor='k', 
                alpha=0.7, 
                s=80)

plt.title('2D PCA of Gene Expression Data')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
plt.grid(True)
plt.tight_layout()
plt.show()

# Correlation heatmap
corr_matrix = final_data.drop(columns=['Cancer_Type']).corr()
sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm', linewidths=0.5)
plt.title("Gene Expression Correlation Heatmap")
plt.show()
