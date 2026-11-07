"""
Predictive Machine Learning Modeling for TCGA-LUAD Patient Survival.

Trains Random Forest classifiers with stratified 5-fold cross validation,
evaluates test ROC-AUC / Confusion Matrix, and computes feature importances.
"""

import os
import argparse
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_validate, train_test_split
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, classification_report, confusion_matrix
)


def prepare_features(df: pd.DataFrame):
    """Preprocess integrated features and extract target labels."""
    df = df.copy()
    if "survival_label" not in df.columns:
        df["survival_label"] = df["Overall Survival Status"].map({
            "0:LIVING": 0,
            "1:DECEASED": 1
        })
    df = df.dropna(subset=["survival_label"])

    exclude = ["SAMPLE_ID", "Overall Survival Status", "survival_label"]
    df_model = df.drop(columns=[c for c in exclude if c in df.columns])

    for col in df_model.select_dtypes(include="object").columns:
        df_model[col] = df_model[col].astype("category").cat.codes

    df_model = df_model.fillna(df_model.median())
    return df_model, df["survival_label"]


def evaluate_model(csv_path: str):
    """Run cross-validation, train/test split, and report metrics."""
    print(f"[*] Loading dataset from: {csv_path}")
    df = pd.read_csv(csv_path)
    X, y = prepare_features(df)

    print(f"    Features shape: {X.shape}, Class distribution: {dict(y.value_counts())}")

    model = RandomForestClassifier(
        n_estimators=200,
        max_depth=8,
        class_weight="balanced",
        random_state=42
    )

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    scores = cross_validate(
        model, X, y, cv=cv,
        scoring=["accuracy", "precision", "recall", "f1", "roc_auc"]
    )

    print("\n--- 5-FOLD STRATIFIED CROSS-VALIDATION RESULTS ---")
    for metric, values in scores.items():
        if "test" in metric:
            name = metric.replace("test_", "").upper()
            print(f"  {name:10s}: {np.mean(values):.4f} (+/- {np.std(values):.4f})")

    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    model.fit(X_train, y_train)
    preds = model.predict(X_test)
    probas = model.predict_proba(X_test)[:, 1]

    print("\n--- HELD-OUT TEST EVALUATION (N = 84) ---")
    print(classification_report(y_test, preds, target_names=["Living (0)", "Deceased (1)"]))
    print(f"Test Accuracy: {accuracy_score(y_test, preds):.4f}")
    print(f"Test Precision:{precision_score(y_test, preds):.4f}")
    print(f"Test Recall:   {recall_score(y_test, preds):.4f}")
    print(f"Test F1-Score: {f1_score(y_test, preds):.4f}")
    print(f"Test ROC-AUC:  {roc_auc_score(y_test, probas):.4f}")

    # Feature Importance
    importances = pd.Series(model.feature_importances_, index=X.columns).sort_values(ascending=False)
    print("\n--- TOP 10 FEATURE IMPORTANCE ---")
    for rank, (feat, score) in enumerate(importances.head(10).items(), 1):
        print(f"  {rank:2d}. {feat:30s}: {score:.4f}")

    return model, importances


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TCGA-LUAD Predictive Modeling")
    parser.add_argument("--data", type=str, default="data/integrated_luad_dataset.csv",
                        help="Path to integrated dataset CSV")
    args = parser.parse_args()

    evaluate_model(args.data)
