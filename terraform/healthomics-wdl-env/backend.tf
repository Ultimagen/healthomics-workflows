terraform {
  backend "s3" {
    # run init.sh to have full backend configuration
    region = "us-east-1"
  }
  required_providers {
    aws = {
      version = ">= 5.43.0"
    }
    local = {
      source = "hashicorp/local"
    }
  }
}